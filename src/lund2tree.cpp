#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Kinematics.h"

const bool fill_original = false;

// A simplified struct holding each final‐state particle’s information
struct SimplePart {
    int id, pid, parentid, parentpid, finalFlag;
    double px, py, pz;
    double E, p, theta, phi;
};

// Macro entry point: convert LUND to multiple dihadron TTrees
void lund2tree(const char* infile, const char* outputFile = "dihadron.root") {
    std::ifstream in(infile);
    if (!in) {
        std::cerr << "Error opening file: " << infile << std::endl;
        return;
    }

    // Create the output ROOT file
    TFile* fOut = new TFile(outputFile, "RECREATE");

    // --- 1) Create the 10 TTrees: original + acceptance for each channel ---
    TTree* tree_ppim        = new TTree("dih_piplus_piminus",            "piplus-piminus dihadron");
    TTree* tree_ppim_acc    = new TTree("dih_piplus_piminus_acceptance", "piplus-piminus dihadron (acceptance)");
    TTree* tree_ppiz        = new TTree("dih_piplus_pi0",                 "piplus-pi0 dihadron");
    TTree* tree_ppiz_acc    = new TTree("dih_piplus_pi0_acceptance",      "piplus-pi0 dihadron (acceptance)");
    TTree* tree_pmiz        = new TTree("dih_piminus_pi0",                "piminus-pi0 dihadron");
    TTree* tree_pmiz_acc    = new TTree("dih_piminus_pi0_acceptance",     "piminus-pi0 dihadron (acceptance)");
    TTree* tree_pppp        = new TTree("dih_piplus_piplus",              "piplus-piplus dihadron");
    TTree* tree_pppp_acc    = new TTree("dih_piplus_piplus_acceptance",   "piplus-piplus dihadron (acceptance)");
    TTree* tree_pimpimp     = new TTree("dih_piminus_piminus",            "piminus-piminus dihadron");
    TTree* tree_pimpimp_acc = new TTree("dih_piminus_piminus_acceptance", "piminus-piminus dihadron (acceptance)");

    // --- 2) Common branch variables for all trees ---
    double x, Q2, y, W, nu;
    int hel;
    double px_e, py_e, pz_e, E_e, p_e, theta_e, phi_e;

    double px_a, py_a, pz_a, E_a, p_a, theta_a, phi_a;
    double px_b, py_b, pz_b, E_b, p_b, theta_b, phi_b;

    double Mh, z1, z2, z, phi_h, phi_R, dihadron_th, pT_lab, Mx, xF1, xF2;
    int is_rho   = false;
    int is_omega = false;
    // Helper lambda to set up identical branches on any TTree pointer
    auto setupBranches = [&](TTree* t) {
        t->Branch("x",   &x,   "x/D");
        t->Branch("Q2",  &Q2,  "Q2/D");
        t->Branch("y",   &y,   "y/D");
        t->Branch("W",   &W,   "W/D");
        t->Branch("nu",  &nu,  "nu/D");
        t->Branch("hel", &hel, "hel/I");

        t->Branch("px_e",   &px_e,   "px_e/D");
        t->Branch("py_e",   &py_e,   "py_e/D");
        t->Branch("pz_e",   &pz_e,   "pz_e/D");
        t->Branch("E_e",    &E_e,    "E_e/D");
        t->Branch("p_e",    &p_e,    "p_e/D");
        t->Branch("theta_e",&theta_e,"theta_e/D");
        t->Branch("phi_e",  &phi_e,  "phi_e/D");

        t->Branch("px_a",    &px_a,    "px_a/D");
        t->Branch("py_a",    &py_a,    "py_a/D");
        t->Branch("pz_a",    &pz_a,    "pz_a/D");
        t->Branch("E_a",     &E_a,     "E_a/D");
        t->Branch("p_a",     &p_a,     "p_a/D");
        t->Branch("theta_a", &theta_a, "theta_a/D");
        t->Branch("phi_a",   &phi_a,   "phi_a/D");

        t->Branch("px_b",    &px_b,    "px_b/D");
        t->Branch("py_b",    &py_b,    "py_b/D");
        t->Branch("pz_b",    &pz_b,    "pz_b/D");
        t->Branch("E_b",     &E_b,     "E_b/D");
        t->Branch("p_b",     &p_b,     "p_b/D");
        t->Branch("theta_b", &theta_b, "theta_b/D");
        t->Branch("phi_b",   &phi_b,   "phi_b/D");

        t->Branch("Mh",      &Mh,         "Mh/D");
        t->Branch("Mx",      &Mx,         "Mx/D");
        t->Branch("xF1",     &xF1,        "xF1/D");
        t->Branch("xF2",     &xF2,        "xF2/D");
        t->Branch("z1",      &z1,         "z1/D");
        t->Branch("z2",      &z2,         "z2/D");
        t->Branch("z",       &z,          "z/D");
        t->Branch("phi_h",   &phi_h,      "phi_h/D");
        t->Branch("phi_R1",  &phi_R,      "phi_R1/D");
        t->Branch("th",      &dihadron_th,"th/D");
        t->Branch("pTtot",  &pT_lab,     "pTtot/D");
        t->Branch("is_rho",   &is_rho,   "is_rho/I");
        t->Branch("is_omega", &is_omega, "is_omega/I");
    };

    // Apply the same branches to all 10 TTrees
    for (auto tptr : { tree_ppim, tree_ppim_acc,
                       tree_ppiz, tree_ppiz_acc,
                       tree_pmiz, tree_pmiz_acc,
                       tree_pppp, tree_pppp_acc,
                       tree_pimpimp, tree_pimpimp_acc })
    {
        setupBranches(tptr);
    }

    // Instantiate your kinematics‐calculator helper
    Kinematics kin;

    // Figure out helicity from filename
    std::string fname(infile);
    if      (fname.find("_LU_p_") != std::string::npos) hel = -1;
    else if (fname.find("_LU_n_") != std::string::npos) hel = +1;
    else                                                    hel =  0;

    // Acceptance angular cuts (in radians)
    const double theta_min = 5 * M_PI / 180.0;
    const double theta_max = 35 * M_PI / 180.0;

    // Loop over each event by reading LUND “header” + nParticles lines
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // Read the LUND header line: nParticles, Mp, targetZ, targetA, hel_ss, beamPid, beamE
        std::istringstream header_ss(line);
        int    nParticles; 
        double Mp; 
        int    targetZ, targetA;
        int    hel_ss, beamPid; 
        double beamE;
        header_ss >> nParticles >> Mp >> targetZ >> targetA >> hel_ss >> beamPid >> beamE;

        // Read the next nParticles lines into a vector<SimplePart>
        std::vector<SimplePart> parts;
        parts.reserve(nParticles);
        for (int i = 0; i < nParticles; ++i) {
            std::getline(in, line);
            if (line.empty()) { --i; continue; }
            std::istringstream iss(line);

            int idx, finalFlag, dummy;
            iss >> idx >> finalFlag >> dummy;
            SimplePart p;
            p.id        = idx;
            p.finalFlag = finalFlag;
            iss >> p.pid;
            iss >> p.parentid >> dummy;
            iss >> p.px >> p.py >> p.pz >> p.E >> dummy;
            p.p = std::hypot(p.px, p.py, p.pz);
            p.theta = (p.p > 0 ? std::acos(p.pz / p.p) : 0);
            p.phi   = std::atan2(p.py, p.px);

            // Find the parent’s PID if parentid > 0
            if (p.parentid > 0) {
                int target_id = p.parentid;
                auto it = std::find_if(parts.begin(), parts.end(),
                                       [target_id](auto const& other){
                                           return other.id == target_id;
                                       });
                p.parentpid = (it != parts.end() ? it->pid : -1);
            }
            else {
                p.parentpid = -1;
            }

            parts.push_back(p);
        } // end of reading all particles

        // --- Identify the highest‐energy final‐state electron (pid = 11, finalFlag==1) ---
        int ie = -1;
        double maxE = -1.0;
        for (int i = 0; i < (int)parts.size(); ++i) {
            if (parts[i].pid == 11 && parts[i].finalFlag == 1 && parts[i].E > maxE) {
                ie = i;
                maxE = parts[i].E;
            }
        }
        if (ie < 0) continue;  // no final‐state electron → skip event

        // Collect charged‐pion indices
        std::vector<int> pip_idx, pim_idx;
        for (int i = 0; i < (int)parts.size(); ++i) {
            if (parts[i].finalFlag != 1) continue;
            if (parts[i].pid == 211)  pip_idx.push_back(i);
            if (parts[i].pid == -211) pim_idx.push_back(i);
        }

        // Collect photon indices whose parentpid == 111 (i.e. π⁰ → γγ)
        std::vector<int> pho_idx;
        for (int i = 0; i < (int)parts.size(); ++i) {
            if (parts[i].pid == 22 && parts[i].parentpid == 111 && parts[i].finalFlag == 1) {
                pho_idx.push_back(i);
            }
        }

        // Precompute four‐vectors for target, beam, and electron
        TLorentzVector lv_target(0, 0, 0, Mp);
        TLorentzVector lv_beam(0, 0, beamE, beamE);
        const auto& e = parts[ie];
        TLorentzVector lv_e(e.px, e.py, e.pz, e.E);

        // A helper lambda that _only_ does the “branch‐variable setting + Fill()”:
        auto fillTree = [&](TTree* tree, const TLorentzVector& lv_a, const TLorentzVector& lv_b){
            // Set the electron branches
            px_e      = e.px;  py_e     = e.py;  pz_e     = e.pz;  E_e  = e.E;
            p_e       = e.p;   theta_e  = e.theta; phi_e   = e.phi;

            // Set particle “a” branches
            px_a      = lv_a.Px(); py_a     = lv_a.Py(); pz_a    = lv_a.Pz();
            E_a       = lv_a.E();  p_a      = lv_a.P();  theta_a = lv_a.Theta(); 
            phi_a     = lv_a.Phi();

            // Set particle “b” branches
            px_b      = lv_b.Px(); py_b     = lv_b.Py(); pz_b    = lv_b.Pz();
            E_b       = lv_b.E();  p_b      = lv_b.P();  theta_b = lv_b.Theta();
            phi_b     = lv_b.Phi();

            // Compute kinematics
            TLorentzVector lv_q     = lv_beam - lv_e;
            TLorentzVector lv_dih   = lv_a + lv_b;

            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q, lv_target );
            W     = kin.W(    Q2, Mp,    nu );
            z1    = kin.z(    lv_target, lv_a,   lv_q );
            z2    = kin.z(    lv_target, lv_b,   lv_q );
            z     = kin.z(    lv_target, lv_dih, lv_q );
            Mh    = lv_dih.M();
            pT_lab = lv_dih.Pt();
            phi_h  = kin.phi_h(  lv_q, lv_e, lv_a, lv_b );
            phi_R  = kin.phi_R(  lv_q, lv_e, lv_a, lv_b );
            dihadron_th = kin.com_th(lv_a, lv_b);
            xF1   = kin.xF( lv_q, lv_a, lv_target, W );
            xF2   = kin.xF( lv_q, lv_b, lv_target, W );
            Mx    = (lv_beam + lv_target - lv_e - lv_a - lv_b).M();

            tree->Fill();
        };

        // ————— π⁺–π⁻ channel —————
        for (int ip : pip_idx) {
            for (int im : pim_idx) {
                const auto& pplus  = parts[ip];
                const auto& pminus = parts[im];
                TLorentzVector lv_pip(pplus.px,  pplus.py,  pplus.pz,  pplus.E);
                TLorentzVector lv_pim(pminus.px, pminus.py, pminus.pz, pminus.E);
                int parent_plus  = parts[ip].parentpid;
                int parent_minus = parts[im].parentpid;
                is_rho   = (parent_plus  == parent_minus && parent_plus  == 113);
                is_omega = (parent_plus  == parent_minus && parent_plus  == 223);
                // 1) Fill the “original” TTree unconditionally
                if(fill_original==true){
                    fillTree(tree_ppim, lv_pip, lv_pim);
                }
                // 2) Now apply “acceptance” cuts (angles + new kinematic cuts)
                if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
                    lv_pip.Theta() >= theta_min && lv_pip.Theta() <= theta_max &&
                    lv_pim.Theta() >= theta_min && lv_pim.Theta() <= theta_max)
                {
                    // Compute all branch variables & kinematics, then test cuts:
                    px_e      = e.px;  py_e     = e.py;  pz_e     = e.pz;  E_e  = e.E;
                    p_e       = e.p;   theta_e  = e.theta; phi_e   = e.phi;

                    px_a      = lv_pip.Px(); py_a     = lv_pip.Py(); pz_a    = lv_pip.Pz();
                    E_a       = lv_pip.E();  p_a      = lv_pip.P();  theta_a = lv_pip.Theta();
                    phi_a     = lv_pip.Phi();

                    px_b      = lv_pim.Px(); py_b     = lv_pim.Py(); pz_b    = lv_pim.Pz();
                    E_b       = lv_pim.E();  p_b      = lv_pim.P();  theta_b = lv_pim.Theta();
                    phi_b     = lv_pim.Phi();

                    TLorentzVector lv_q    = lv_beam - lv_e;
                    TLorentzVector lv_dih  = lv_pip + lv_pim;

                    Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
                    y     = kin.y(  beamE, E_e );
                    nu    = kin.nu( beamE, E_e );
                    x     = kin.x(    Q2, lv_q, lv_target );
                    W     = kin.W(    Q2, Mp,    nu );
                    z1    = kin.z(    lv_target, lv_pip, lv_q );
                    z2    = kin.z(    lv_target, lv_pim, lv_q );
                    z     = kin.z(    lv_target, lv_dih, lv_q );
                    Mh    = lv_dih.M();
                    pT_lab = lv_dih.Pt();
                    phi_h  = kin.phi_h(  lv_q, lv_e, lv_pip, lv_pim );
                    phi_R  = kin.phi_R(  lv_q, lv_e, lv_pip, lv_pim );
                    dihadron_th = kin.com_th(lv_pip, lv_pim);
                    xF1   = kin.xF( lv_q, lv_pip, lv_target, W );
                    xF2   = kin.xF( lv_q, lv_pim, lv_target, W );
                    Mx    = (lv_beam + lv_target - lv_e - lv_pip - lv_pim).M();

                    // Apply ALL extra cuts for π⁺–π⁻ (including p_b > 1.25):
                    if (Mx  > 1.5 &&
                        xF1 > 0   &&
                        xF2 > 0   &&
                        z   < 0.95&&
                        W   > 2   &&
                        Q2  > 1   &&
                        p_a > 1.25&&
                        p_b > 1.25)
                    {
                        tree_ppim_acc->Fill();
                    }
                }
            }
        }

        // ————— π⁺–π⁰ channel —————
        for (int ip : pip_idx) {
            const auto& pplus = parts[ip];
            TLorentzVector lv_pip(pplus.px, pplus.py, pplus.pz, pplus.E);

            // Loop over all distinct γγ pairs reconstructing π⁰
            for (size_t i1 = 0; i1 < pho_idx.size(); ++i1) {
                for (size_t i2 = i1 + 1; i2 < pho_idx.size(); ++i2) {
                    int i = pho_idx[i1], j = pho_idx[i2];
                    if (parts[i].parentid != parts[j].parentid) continue; // require same parent (=111)
                    
                    const auto& ph1 = parts[i];
                    const auto& ph2 = parts[j];

                    // Enforce each photon’s energy ≥ 0.2 GeV
                    if (ph1.E < 0.2 || ph2.E < 0.2) continue;

                    TLorentzVector lv_ph1(ph1.px, ph1.py, ph1.pz, ph1.E);
                    TLorentzVector lv_ph2(ph2.px, ph2.py, ph2.pz, ph2.E);
                    TLorentzVector lv_pi0 = lv_ph1 + lv_ph2;
                    int pi0_id = parts[i].parentid;  
                    auto it0   = std::find_if(parts.begin(), parts.end(),
                                  [pi0_id](auto const& p){ return p.id == pi0_id; });
                    int parent_pi0 = (it0 != parts.end() ? it0->parentpid : -1);
                    int parent_pip = parts[ip].parentpid;
                    is_rho   = (parent_pip == 213 && parent_pi0 == 213);
                    is_omega = (parent_pip == 223 && parent_pi0 == 223);
                    // 1) Fill the “original” π⁺–π⁰ tree (no acceptance cuts)
                    if(fill_original==true){
                        fillTree(tree_ppiz, lv_pip, lv_pi0);
                    }
                    // 2) Now acceptance angles + kinematic cuts (but **no** p_b > 1.25):
                    if (lv_e.Theta()  >= theta_min && lv_e.Theta()  <= theta_max &&
                        lv_pip.Theta() >= theta_min && lv_pip.Theta() <= theta_max &&
                        lv_ph1.Theta()  >= theta_min && lv_ph1.Theta()  <= theta_max &&
                        lv_ph2.Theta()  >= theta_min && lv_ph2.Theta()  <= theta_max)
                    {
                        // Compute branch variables + kinematics
                        px_e      = e.px;  py_e     = e.py;  pz_e     = e.pz;  E_e  = e.E;
                        p_e       = e.p;   theta_e  = e.theta; phi_e   = e.phi;

                        px_a      = lv_pip.Px(); py_a  = lv_pip.Py(); pz_a   = lv_pip.Pz();
                        E_a       = lv_pip.E();  p_a   = lv_pip.P();  theta_a= lv_pip.Theta();
                        phi_a     = lv_pip.Phi();

                        px_b      = lv_pi0.Px(); py_b  = lv_pi0.Py(); pz_b   = lv_pi0.Pz();
                        E_b       = lv_pi0.E();  p_b   = lv_pi0.P();  theta_b= lv_pi0.Theta();
                        phi_b     = lv_pi0.Phi();

                        TLorentzVector lv_q    = lv_beam - lv_e;
                        TLorentzVector lv_dih  = lv_pip + lv_pi0;

                        Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
                        y     = kin.y(  beamE, E_e );
                        nu    = kin.nu( beamE, E_e );
                        x     = kin.x(    Q2, lv_q, lv_target );
                        W     = kin.W(    Q2, Mp,    nu );
                        z1    = kin.z(    lv_target, lv_pip, lv_q );
                        z2    = kin.z(    lv_target, lv_pi0, lv_q );
                        z     = kin.z(    lv_target, lv_dih, lv_q );
                        Mh    = lv_dih.M();
                        pT_lab = lv_dih.Pt();
                        phi_h  = kin.phi_h(  lv_q, lv_e, lv_pip, lv_pi0 );
                        phi_R  = kin.phi_R(  lv_q, lv_e, lv_pip, lv_pi0 );
                        dihadron_th = kin.com_th(lv_pip, lv_pi0);
                        xF1   = kin.xF( lv_q, lv_pip,  lv_target, W );
                        xF2   = kin.xF( lv_q, lv_pi0,  lv_target, W );
                        Mx    = (lv_beam + lv_target - lv_e - lv_pip - lv_pi0).M();

                        // Now apply π⁺–π⁰ acceptance cuts **excluding** p_b
                        if (Mx  > 1.5 &&
                            xF1 > 0   &&
                            xF2 > 0   &&
                            z   < 0.95&&
                            W   > 2   &&
                            Q2  > 1   &&
                            p_a > 1.25)
                        {
                            tree_ppiz_acc->Fill();
                        }
                    }
                }
            }
        }

        // ————— π⁻–π⁰ channel —————
        for (int im : pim_idx) {
            const auto& pminus = parts[im];
            TLorentzVector lv_pim(pminus.px, pminus.py, pminus.pz, pminus.E);

            for (size_t i1 = 0; i1 < pho_idx.size(); ++i1) {
                for (size_t i2 = i1 + 1; i2 < pho_idx.size(); ++i2) {
                    int i = pho_idx[i1], j = pho_idx[i2];
                    if (parts[i].parentid != parts[j].parentid) continue;

                    const auto& ph1 = parts[i];
                    const auto& ph2 = parts[j];

                    // Require each photon energy ≥ 0.2 GeV
                    if (ph1.E < 0.2 || ph2.E < 0.2) continue;

                    TLorentzVector lv_ph1(ph1.px, ph1.py, ph1.pz, ph1.E);
                    TLorentzVector lv_ph2(ph2.px, ph2.py, ph2.pz, ph2.E);
                    TLorentzVector lv_pi0 = lv_ph1 + lv_ph2;
                    int pi0_id = parts[i].parentid;  
                    auto it0   = std::find_if(parts.begin(), parts.end(),
                                  [pi0_id](auto const& p){ return p.id == pi0_id; });
                    int parent_pi0 = (it0 != parts.end() ? it0->parentpid : -1);
                    int parent_pim = parts[im].parentpid;
                    is_rho   = (parent_pim == -213 && parent_pi0 == -213);
                    is_omega = (parent_pim ==  223 && parent_pi0 ==  223);
                    // Fill original π⁻–π⁰ tree
                    if(fill_original==true){
                        fillTree(tree_pmiz, lv_pim, lv_pi0);
                    }
                    // Acceptance angles + kinematic cuts (no p_b ):
                    if (lv_e.Theta()  >= theta_min && lv_e.Theta()  <= theta_max &&
                        lv_pim.Theta() >= theta_min && lv_pim.Theta() <= theta_max &&
                        lv_ph1.Theta()  >= theta_min && lv_ph1.Theta()  <= theta_max &&
                        lv_ph2.Theta()  >= theta_min && lv_ph2.Theta()  <= theta_max)
                    {
                        // Compute all branches & kinematics
                        px_e      = e.px;  py_e     = e.py;  pz_e     = e.pz;  E_e  = e.E;
                        p_e       = e.p;   theta_e  = e.theta; phi_e   = e.phi;

                        px_a      = lv_pim.Px(); py_a   = lv_pim.Py(); pz_a   = lv_pim.Pz();
                        E_a       = lv_pim.E();  p_a    = lv_pim.P();  theta_a= lv_pim.Theta();
                        phi_a     = lv_pim.Phi();

                        px_b      = lv_pi0.Px(); py_b   = lv_pi0.Py(); pz_b   = lv_pi0.Pz();
                        E_b       = lv_pi0.E();  p_b    = lv_pi0.P();  theta_b= lv_pi0.Theta();
                        phi_b     = lv_pi0.Phi();

                        TLorentzVector lv_q    = lv_beam - lv_e;
                        TLorentzVector lv_dih  = lv_pim + lv_pi0;

                        Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
                        y     = kin.y(  beamE, E_e );
                        nu    = kin.nu( beamE, E_e );
                        x     = kin.x(    Q2, lv_q, lv_target );
                        W     = kin.W(    Q2, Mp,    nu );
                        z1    = kin.z(    lv_target, lv_pim, lv_q );
                        z2    = kin.z(    lv_target, lv_pi0, lv_q );
                        z     = kin.z(    lv_target, lv_dih, lv_q );
                        Mh    = lv_dih.M();
                        pT_lab = lv_dih.Pt();
                        phi_h  = kin.phi_h(  lv_q, lv_e, lv_pim, lv_pi0 );
                        phi_R  = kin.phi_R(  lv_q, lv_e, lv_pim, lv_pi0 );
                        dihadron_th = kin.com_th(lv_pim, lv_pi0);
                        xF1   = kin.xF( lv_q, lv_pim,  lv_target, W );
                        xF2   = kin.xF( lv_q, lv_pi0,  lv_target, W );
                        Mx    = (lv_beam + lv_target - lv_e - lv_pim - lv_pi0).M();

                        // Apply π⁻–π⁰ acceptance cuts (excluding p_b):
                        if (Mx  > 1.5 &&
                            xF1 > 0   &&
                            xF2 > 0   &&
                            z   < 0.95&&
                            W   > 2   &&
                            Q2  > 1   &&
                            p_a > 1.25)
                        {
                            tree_pmiz_acc->Fill();
                        }
                    }
                }
            }
        }

        // ————— π⁺–π⁺ channel —————
        for (size_t i1 = 0; i1 < pip_idx.size(); ++i1) {
            for (size_t i2 = i1 + 1; i2 < pip_idx.size(); ++i2) {
                int ia = pip_idx[i1], ib = pip_idx[i2];
                TLorentzVector lv1(parts[ia].px, parts[ia].py, parts[ia].pz, parts[ia].E);
                TLorentzVector lv2(parts[ib].px, parts[ib].py, parts[ib].pz, parts[ib].E);

                // Determine leading (lo) / subleading (hi) by z:
                TLorentzVector lv_q   = lv_beam - lv_e;
                double z1_ = kin.z(lv_target, lv1, lv_q);
                double z2_ = kin.z(lv_target, lv2, lv_q);
                TLorentzVector lo = (z1_ >= z2_ ? lv1 : lv2);
                TLorentzVector hi = (z1_ >= z2_ ? lv2 : lv1);

                // 1) Fill the “original” π⁺–π⁺ tree
                if(fill_original==true){
                    fillTree(tree_pppp, lo, hi);
                }
                // 2) Acceptance + cuts (including p_b > 1.25):
                if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
                    lo.Theta()   >= theta_min && lo.Theta()   <= theta_max &&
                    hi.Theta()   >= theta_min && hi.Theta()   <= theta_max)
                {
                    // Compute all branches & kinematics
                    px_e      = e.px;  py_e     = e.py;  pz_e     = e.pz;  E_e  = e.E;
                    p_e       = e.p;   theta_e  = e.theta; phi_e   = e.phi;

                    px_a      = lo.Px(); py_a   = lo.Py(); pz_a   = lo.Pz();
                    E_a       = lo.E();  p_a    = lo.P();  theta_a= lo.Theta();
                    phi_a     = lo.Phi();

                    px_b      = hi.Px(); py_b   = hi.Py(); pz_b   = hi.Pz();
                    E_b       = hi.E();  p_b    = hi.P();  theta_b= hi.Theta();
                    phi_b     = hi.Phi();

                    TLorentzVector lv_q2   = lv_beam - lv_e;
                    TLorentzVector lv_dih2 = lo + hi;

                    Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
                    y     = kin.y(  beamE, E_e );
                    nu    = kin.nu( beamE, E_e );
                    x     = kin.x(    Q2, lv_q2, lv_target );
                    W     = kin.W(    Q2, Mp,    nu );
                    z1    = kin.z(    lv_target, lo,    lv_q2 );
                    z2    = kin.z(    lv_target, hi,    lv_q2 );
                    z     = kin.z(    lv_target, lv_dih2, lv_q2 );
                    Mh    = lv_dih2.M();
                    pT_lab = lv_dih2.Pt();
                    phi_h  = kin.phi_h(  lv_q2, lv_e, lo, hi );
                    phi_R  = kin.phi_R(  lv_q2, lv_e, lo, hi );
                    dihadron_th = kin.com_th(lo, hi);
                    xF1   = kin.xF( lv_q2, lo,  lv_target, W );
                    xF2   = kin.xF( lv_q2, hi,  lv_target, W );
                    Mx    = (lv_beam + lv_target - lv_e - lo - hi).M();

                    if (Mx  > 1.5 &&
                        xF1 > 0   &&
                        xF2 > 0   &&
                        z   < 0.95&&
                        W   > 2   &&
                        Q2  > 1   &&
                        p_a > 1.25&&
                        p_b > 1.25)
                    {
                        tree_pppp_acc->Fill();
                    }
                }
            }
        }

        // ————— π⁻–π⁻ channel —————
        for (size_t i1 = 0; i1 < pim_idx.size(); ++i1) {
            for (size_t i2 = i1 + 1; i2 < pim_idx.size(); ++i2) {
                int ia = pim_idx[i1], ib = pim_idx[i2];
                TLorentzVector lv1(parts[ia].px, parts[ia].py, parts[ia].pz, parts[ia].E);
                TLorentzVector lv2(parts[ib].px, parts[ib].py, parts[ib].pz, parts[ib].E);

                TLorentzVector lv_q2   = lv_beam - lv_e;
                double z1_ = kin.z(lv_target, lv1, lv_q2);
                double z2_ = kin.z(lv_target, lv2, lv_q2);
                TLorentzVector lo = (z1_ >= z2_ ? lv1 : lv2);
                TLorentzVector hi = (z1_ >= z2_ ? lv2 : lv1);

                // 1) Fill original π⁻–π⁻ tree
                if(fill_original==true){
                    fillTree(tree_pimpimp, lo, hi);
                }
                // 2) Acceptance + cuts (including p_b > 1.25):
                if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
                    lo.Theta()   >= theta_min && lo.Theta()   <= theta_max &&
                    hi.Theta()   >= theta_min && hi.Theta()   <= theta_max)
                {
                    px_e      = e.px;  py_e     = e.py;  pz_e     = e.pz;  E_e  = e.E;
                    p_e       = e.p;   theta_e  = e.theta; phi_e   = e.phi;

                    px_a      = lo.Px(); py_a   = lo.Py(); pz_a   = lo.Pz();
                    E_a       = lo.E();  p_a    = lo.P();  theta_a= lo.Theta();
                    phi_a     = lo.Phi();

                    px_b      = hi.Px(); py_b   = hi.Py(); pz_b   = hi.Pz();
                    E_b       = hi.E();  p_b    = hi.P();  theta_b= hi.Theta();
                    phi_b     = hi.Phi();

                    TLorentzVector lv_q3   = lv_beam - lv_e;
                    TLorentzVector lv_dih3 = lo + hi;

                    Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
                    y     = kin.y(  beamE, E_e );
                    nu    = kin.nu( beamE, E_e );
                    x     = kin.x(    Q2, lv_q3, lv_target );
                    W     = kin.W(    Q2, Mp,    nu );
                    z1    = kin.z(    lv_target, lo,     lv_q3 );
                    z2    = kin.z(    lv_target, hi,     lv_q3 );
                    z     = kin.z(    lv_target, lv_dih3, lv_q3 );
                    Mh    = lv_dih3.M();
                    pT_lab = lv_dih3.Pt();
                    phi_h  = kin.phi_h(  lv_q3, lv_e, lo, hi );
                    phi_R  = kin.phi_R(  lv_q3, lv_e, lo, hi );
                    dihadron_th = kin.com_th(lo, hi);
                    xF1   = kin.xF( lv_q3, lo,  lv_target, W );
                    xF2   = kin.xF( lv_q3, hi,  lv_target, W );
                    Mx    = (lv_beam + lv_target - lv_e - lo - hi).M();

                    if (Mx  > 1.5 &&
                        xF1 > 0   &&
                        xF2 > 0   &&
                        z   < 0.95&&
                        W   > 2   &&
                        Q2  > 1   &&
                        p_a > 1.25&&
                        p_b > 1.25)
                    {
                        tree_pimpimp_acc->Fill();
                    }
                }
            }
        }

    } // end of event loop

    // Write all TTrees to disk
    fOut->Write();
    fOut->Close();
}

// --------------------------- main() ---------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.dat> [output.root]\n";
        return 1;
    }
    lund2tree(argv[1], (argc > 2 ? argv[2] : "dihadron.root"));
    return 0;
}







// #include <TFile.h>
// #include <TTree.h>
// #include <TLorentzVector.h>
// #include <fstream>
// #include <sstream>
// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <algorithm>
// #include "Kinematics.h"

// struct SimplePart {
//     int id, pid, parentid, parentpid, finalFlag;
//     double px, py, pz;
//     double E, p, theta, phi;
// };

// // Macro entry point: convert LUND to multiple dihadron TTrees
// void lund2tree(const char* infile, const char* outputFile = "dihadron.root") {
//     std::ifstream in(infile);
//     if (!in) { std::cerr << "Error opening file: " << infile << std::endl; return; }

//     TFile* fOut = new TFile(outputFile, "RECREATE");
    
//     // --- 1) Create the 10 TTrees: original + acceptance for each channel ---
//     TTree* tree_ppim        = new TTree("dih_piplus_piminus",           "piplus-piminus dihadron");
//     TTree* tree_ppim_acc    = new TTree("dih_piplus_piminus_acceptance","piplus-piminus dihadron (acceptance)");
//     TTree* tree_ppiz        = new TTree("dih_piplus_pi0",                "piplus-pi0 dihadron");
//     TTree* tree_ppiz_acc    = new TTree("dih_piplus_pi0_acceptance",     "piplus-pi0 dihadron (acceptance)");
//     TTree* tree_pmiz        = new TTree("dih_piminus_pi0",               "piminus-pi0 dihadron");
//     TTree* tree_pmiz_acc    = new TTree("dih_piminus_pi0_acceptance",    "piminus-pi0 dihadron (acceptance)");
//     TTree* tree_pppp        = new TTree("dih_piplus_piplus",             "piplus-piplus dihadron");
//     TTree* tree_pppp_acc    = new TTree("dih_piplus_piplus_acceptance",  "piplus-piplus dihadron (acceptance)");
//     TTree* tree_pimpimp     = new TTree("dih_piminus_piminus",           "piminus-piminus dihadron");
//     TTree* tree_pimpimp_acc = new TTree("dih_piminus_piminus_acceptance","piminus-piminus dihadron (acceptance)");

//     // Common branches
//     double x, Q2, y, W, nu;
//     int hel;
//     double px_e, py_e, pz_e, E_e, p_e, theta_e, phi_e;

//     double px_a, py_a, pz_a, E_a, p_a, theta_a, phi_a;
//     double px_b, py_b, pz_b, E_b, p_b, theta_b, phi_b;

//     double Mh, z1, z2, z, phi_h, phi_R, dihadron_th, pT_lab, Mx, xF1, xF2;

//     // Set branches on each tree
//     auto setupBranches = [&](TTree* t) {
//         t->Branch("x", &x, "x/D");
//         t->Branch("Q2", &Q2, "Q2/D");
//         t->Branch("y", &y, "y/D");
//         t->Branch("W", &W, "W/D");
//         t->Branch("nu", &nu, "nu/D");
//         t->Branch("hel", &hel, "hel/I");
//         t->Branch("px_e", &px_e, "px_e/D");
//         t->Branch("py_e", &py_e, "py_e/D");
//         t->Branch("pz_e", &pz_e, "pz_e/D");
//         t->Branch("E_e",  &E_e,  "E_e/D");
//         t->Branch("p_e",  &p_e,  "p_e/D");
//         t->Branch("theta_e", &theta_e, "theta_e/D");
//         t->Branch("phi_e",   &phi_e,   "phi_e/D");

//         t->Branch("px_a",  &px_a,  "px_a/D");
//         t->Branch("py_a",  &py_a,  "py_a/D");
//         t->Branch("pz_a",  &pz_b,  "pz_a/D");
//         t->Branch("E_a",   &E_a,   "E_a/D");
//         t->Branch("p_a",   &p_a,   "p_a/D");
//         t->Branch("theta_a", &theta_a, "theta_a/D");
//         t->Branch("phi_a",   &phi_a,   "phi_a/D");

//         t->Branch("px_b",  &px_b,  "px_b/D");
//         t->Branch("py_b",  &py_b,  "py_b/D");
//         t->Branch("pz_b",  &pz_b,  "pz_b/D");
//         t->Branch("E_b",   &E_b,   "E_b/D");
//         t->Branch("p_b",   &p_b,   "p_b/D");
//         t->Branch("theta_b", &theta_b, "theta_b/D");
//         t->Branch("phi_b",   &phi_b,   "phi_b/D");

//         t->Branch("Mh",      &Mh, "Mh/D");
//         t->Branch("Mx",      &Mx, "Mx/D");
//         t->Branch("xF1",      &xF1, "xF1/D");
//         t->Branch("xF2",      &xF2, "xF2/D");
//         t->Branch("z1",      &z1, "z1/D");
//         t->Branch("z2",      &z2, "z2/D");
//         t->Branch("z",       &z,  "z/D");
//         t->Branch("phi_h",   &phi_h, "phi_h/D");
//         t->Branch("phi_R1",   &phi_R, "phi_R1/D");
//         t->Branch("th",      &dihadron_th, "th/D");
//         t->Branch("pT_lab",  &pT_lab,       "pT_lab/D");
//     };

//     // apply to all trees
//     for (auto t : { tree_ppim, tree_ppim_acc,
//                     tree_ppiz, tree_ppiz_acc,
//                     tree_pmiz, tree_pmiz_acc,
//                     tree_pppp, tree_pppp_acc,
//                     tree_pimpimp, tree_pimpimp_acc })
//     {
//         setupBranches(t);
//     }
    
//     Kinematics kin;
//     std::string line;
//     std::string fname(infile);
//     if      (fname.find("_LU_p_")!=std::string::npos) hel = -1; 
//     else if (fname.find("_LU_n_")!=std::string::npos) hel = +1;
//     else                                                 hel =  0;

//     // acceptance angles in radians
//     const double theta_min = 5 * M_PI / 180.0;
//     const double theta_max = 35 * M_PI / 180.0;
    
//     while (std::getline(in, line)) {
//         if (line.empty()) continue;
//         std::istringstream header_ss(line);
//         int nParticles; double Mp; int targetZ,targetA; int hel_ss, beamPid; double beamE;
//         header_ss >> nParticles >> Mp >> targetZ >> targetA >> hel_ss >> beamPid >> beamE;
        
//         std::vector<SimplePart> parts;
//         parts.reserve(nParticles);
//         for (int i=0; i<nParticles; ++i) {
//             std::getline(in, line);
//             if (line.empty()) { --i; continue; }
//             std::istringstream iss(line);
//             int idx, finalFlag, dummy; iss >> idx >> finalFlag >> dummy;
//             // if (finalFlag!=1) continue;
//             SimplePart p;
//             p.finalFlag = finalFlag;
//             p.id = idx;
//             iss >> p.pid;
//             iss >> p.parentid >> dummy;
//             iss >> p.px >> p.py >> p.pz >> p.E >> dummy;
//             p.p = std::hypot(p.px,p.py,p.pz);
//             p.theta = (p.p>0 ? std::acos(p.pz/p.p) : 0);
//             p.phi = std::atan2(p.py,p.px);
//             if (p.parentid > 0) {
//                 int target_id = p.parentid;
//                 auto it = std::find_if(parts.begin(), parts.end(),
//                     [target_id](auto const& other){
//                         return other.id == target_id;
//                     });
//                 if (it != parts.end()) {
//                     p.parentpid = it->pid;
//                 }
//                 else {
//                     // no matching mother found
//                     p.parentpid = -1;
//                 }
//             }
//             else {
//                 p.parentpid = -1;
//             }
//             parts.push_back(p);
//         }
//         // find highest-energy electron
//         int ie=-1; double maxE=-1;
//         for (int i=0; i<parts.size(); ++i) if (parts[i].pid==11 && parts[i].E>maxE && parts[i].finalFlag==1) { ie=i; maxE=parts[i].E; }
//         if (ie<0) continue;

//         // collect charged pions
//         std::vector<int> pip_idx, pim_idx;
//         for (int i=0; i<parts.size(); ++i) {
//             if (parts[i].finalFlag!=1) continue;
//             if (parts[i].pid==211)  pip_idx.push_back(i);
//             if (parts[i].pid==-211) pim_idx.push_back(i);
//         }
        
//         // collect photons for pi0
//         std::vector<int> pho_idx;
//         for (int i=0; i<parts.size(); ++i) if (parts[i].pid==22 && parts[i].parentpid==111 && parts[i].finalFlag==1) pho_idx.push_back(i);
        
//         // precompute electron lv
//         TLorentzVector lv_target(0,0,0,Mp);
//         TLorentzVector lv_beam(0,0,beamE,beamE);
//         const auto& e = parts[ie];
//         TLorentzVector lv_e(e.px,e.py,e.pz,e.E);
        
//         // helper for kinematics and filling
//         auto fill_and_compute = [&](TTree* tree, const TLorentzVector& lv_a, const TLorentzVector& lv_b) {
//             // electron
//             px_e=e.px; py_e=e.py; pz_e=e.pz; E_e=e.E; p_e=e.p; theta_e=e.theta; phi_e=e.phi;
//             // partner a
//             px_a=lv_a.Px(); py_a=lv_a.Py(); pz_a=lv_a.Pz(); E_a=lv_a.E(); p_a=lv_a.P(); theta_a=lv_a.Theta(); phi_a=lv_a.Phi();
//             // partner b
//             px_b=lv_b.Px(); py_b=lv_b.Py(); pz_b=lv_b.Pz(); E_b=lv_b.E(); p_b=lv_b.P(); theta_b=lv_b.Theta(); phi_b=lv_b.Phi();
            
//             TLorentzVector lv_q = lv_beam - lv_e;
//             TLorentzVector lv_dih = lv_a + lv_b;
//             Q2 = kin.Q2(beamE, E_e, std::cos(theta_e));
//             y  = kin.y(beamE, E_e);
//             nu = kin.nu(beamE, E_e);
//             x  = kin.x(    Q2, lv_q,    lv_target);
//             W  = kin.W(    Q2, Mp,      nu);
//             z1 = kin.z(    lv_target, lv_a, lv_q);
//             z2 = kin.z(    lv_target, lv_b, lv_q);
//             z  = kin.z(    lv_target, lv_dih, lv_q);
//             Mh = lv_dih.M();
//             pT_lab = lv_dih.Pt();
//             phi_h  = kin.phi_h(  lv_q, lv_e, lv_a, lv_b);
//             phi_R  = kin.phi_R(  lv_q, lv_e, lv_a, lv_b);
//             dihadron_th = kin.com_th(lv_a, lv_b);
//             xF1 = kin.xF( lv_q, lv_a, lv_target, W);
//             xF2 = kin.xF( lv_q, lv_b, lv_target, W);
//             Mx = (lv_beam + lv_target - lv_e - lv_a - lv_b).M();
//             tree->Fill();
//         };
        
//         // ---- piplus-piminus ----
//         for (int ip : pip_idx) for (int im : pim_idx) {
//             TLorentzVector lv_pip(parts[ip].px, parts[ip].py, parts[ip].pz, parts[ip].E);
//             TLorentzVector lv_pim(parts[im].px, parts[im].py, parts[im].pz, parts[im].E);

//             // fill original
//             fill_and_compute(tree_ppim, lv_pip, lv_pim);

//             // acceptance cut on electron, pip, pim
//             if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
//                 lv_pip.Theta() >= theta_min && lv_pip.Theta() <= theta_max &&
//                 lv_pim.Theta() >= theta_min && lv_pim.Theta() <= theta_max)
//             {
//                 fill_and_compute(tree_ppim_acc, lv_pip, lv_pim);
//             }
//         }

//         // ---- piplus-pi0 ----
//         for (int ip : pip_idx) {
//             TLorentzVector lv_pip(parts[ip].px, parts[ip].py, parts[ip].pz, parts[ip].E);
//             for (size_t i1 = 0; i1 < pho_idx.size(); ++i1) {
//                 for (size_t i2 = i1 + 1; i2 < pho_idx.size(); ++i2) {
//                     int i = pho_idx[i1], j = pho_idx[i2];
//                     if (parts[i].parentid != parts[j].parentid) continue;

//                     TLorentzVector lv_ph1(parts[i].px, parts[i].py, parts[i].pz, parts[i].E);
//                     TLorentzVector lv_ph2(parts[j].px, parts[j].py, parts[j].pz, parts[j].E);
//                     TLorentzVector lv_pi0 = lv_ph1 + lv_ph2;

//                     // fill original
//                     fill_and_compute(tree_ppiz, lv_pip, lv_pi0);

//                     // acceptance: electron, pip, both photons
//                     if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
//                         lv_pip.Theta() >= theta_min && lv_pip.Theta() <= theta_max &&
//                         lv_ph1.Theta() >= theta_min && lv_ph1.Theta() <= theta_max &&
//                         lv_ph2.Theta() >= theta_min && lv_ph2.Theta() <= theta_max)
//                     {
//                         fill_and_compute(tree_ppiz_acc, lv_pip, lv_pi0);
//                     }
//                 }
//             }
//         }

//         // ---- piminus-pi0 ----
//         for (int im : pim_idx) {
//             TLorentzVector lv_pim(parts[im].px, parts[im].py, parts[im].pz, parts[im].E);
//             for (size_t i1 = 0; i1 < pho_idx.size(); ++i1) {
//                 for (size_t i2 = i1 + 1; i2 < pho_idx.size(); ++i2) {
//                     int i = pho_idx[i1], j = pho_idx[i2];
//                     if (parts[i].parentid != parts[j].parentid) continue;

//                     TLorentzVector lv_ph1(parts[i].px, parts[i].py, parts[i].pz, parts[i].E);
//                     TLorentzVector lv_ph2(parts[j].px, parts[j].py, parts[j].pz, parts[j].E);
//                     TLorentzVector lv_pi0 = lv_ph1 + lv_ph2;

//                     fill_and_compute(tree_pmiz, lv_pim, lv_pi0);

//                     if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
//                         lv_pim.Theta() >= theta_min && lv_pim.Theta() <= theta_max &&
//                         lv_ph1.Theta() >= theta_min && lv_ph1.Theta() <= theta_max &&
//                         lv_ph2.Theta() >= theta_min && lv_ph2.Theta() <= theta_max)
//                     {
//                         fill_and_compute(tree_pmiz_acc, lv_pim, lv_pi0);
//                     }
//                 }
//             }
//         }

//         // ---- piplus-piplus ----
//         for (size_t i1 = 0; i1 < pip_idx.size(); ++i1) {
//             for (size_t i2 = i1 + 1; i2 < pip_idx.size(); ++i2) {
//                 int ia = pip_idx[i1], ib = pip_idx[i2];
//                 TLorentzVector lv1(parts[ia].px, parts[ia].py, parts[ia].pz, parts[ia].E);
//                 TLorentzVector lv2(parts[ib].px, parts[ib].py, parts[ib].pz, parts[ib].E);

//                 // determine leading by z
//                 TLorentzVector lv_q = lv_beam - lv_e;
//                 double z1_ = kin.z(lv_target, lv1, lv_q);
//                 double z2_ = kin.z(lv_target, lv2, lv_q);
//                 auto [lo, hi] = (z1_ >= z2_)
//                     ? std::pair{lv1, lv2}
//                     : std::pair{lv2, lv1};

//                 fill_and_compute(tree_pppp, lo, hi);

//                 if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
//                     lo.Theta()  >= theta_min && lo.Theta()  <= theta_max &&
//                     hi.Theta()  >= theta_min && hi.Theta()  <= theta_max)
//                 {
//                     fill_and_compute(tree_pppp_acc, lo, hi);
//                 }
//             }
//         }

//         // ---- piminus-piminus ----
//         for (size_t i1 = 0; i1 < pim_idx.size(); ++i1) {
//             for (size_t i2 = i1 + 1; i2 < pim_idx.size(); ++i2) {
//                 int ia = pim_idx[i1], ib = pim_idx[i2];
//                 TLorentzVector lv1(parts[ia].px, parts[ia].py, parts[ia].pz, parts[ia].E);
//                 TLorentzVector lv2(parts[ib].px, parts[ib].py, parts[ib].pz, parts[ib].E);

//                 TLorentzVector lv_q = lv_beam - lv_e;
//                 double z1_ = kin.z(lv_target, lv1, lv_q);
//                 double z2_ = kin.z(lv_target, lv2, lv_q);
//                 auto [lo, hi] = (z1_ >= z2_)
//                     ? std::pair{lv1, lv2}
//                     : std::pair{lv2, lv1};

//                 fill_and_compute(tree_pimpimp, lo, hi);

//                 if (lv_e.Theta() >= theta_min && lv_e.Theta() <= theta_max &&
//                     lo.Theta()  >= theta_min && lo.Theta()  <= theta_max &&
//                     hi.Theta()  >= theta_min && hi.Theta()  <= theta_max)
//                 {
//                     fill_and_compute(tree_pimpimp_acc, lo, hi);
//                 }
//             }
//         }
//     }

//     fOut->Write();
//     fOut->Close();
// }

// int main(int argc, char** argv) {
//     if (argc < 2) { std::cerr << "Usage: "<<argv[0]<<" <input.dat> [output.root]\n"; return 1; }
//     lund2tree(argv[1], (argc>2?argv[2]:"dihadron.root"));
//     return 0;
// }
