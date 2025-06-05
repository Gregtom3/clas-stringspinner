#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Pythia8/Event.h>
#include "Tools.h"
#include "EventObjects.h"

namespace clas {

  /// @brief checklist for particle cuts
  class CheckList {

    public:

      /// modes for `CheckList` operations
      enum Mode {
        /// min and max values are not used
        kNoCuts,
        /// cuts applied to list of single hadrons
        k1hCuts,
        /// cuts applied to pair of hadrons
        k2hCuts
      };

      /// @param opt_name option name from command line
      /// @param the mode for this checklist
      CheckList(std::string const& opt_name, Mode mode) :
        m_opt_name(opt_name),
        m_mode(mode),
        m_enabled(false) {}

      /// @brief setup the checklist by parsing command line arguments
      /// @param opt_arg argument string from command line
      void Setup(char const* opt_arg)
      {
        // parse command line arguments
        switch(m_mode) {
          case kNoCuts:
            {
              Tokenize(opt_arg, [&](auto token, auto i) {
                m_pdg_list.push_back(std::stoi(token));
              });
              if(m_pdg_list.empty())
                throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 1 argument", m_opt_name));
              break;
            }
          case k1hCuts:
          case k2hCuts:
            {
              Tokenize(opt_arg, [&](auto token, auto i) {
                switch(i) {
                  case 0: m_min = std::stod(token); break;
                  case 1: m_max = std::stod(token); break;
                  default: m_pdg_list.push_back(std::stoi(token));
                }
              });
              if(m_mode == k2hCuts && m_pdg_list.size() != 2)
                throw std::runtime_error(fmt::format("value of option '--{}' must have exactly 2 PDG values, which specify a dihadron", m_opt_name));
              if(m_pdg_list.empty())
                throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 3 arguments", m_opt_name));
              if(m_min >= m_max)
                throw std::runtime_error(fmt::format("option '--{}' has MIN >= MAX", m_opt_name));
              break;
            }
        }
        m_enabled = true;
      }

      // @returns true if this `CheckList` is enabled
      bool const& Enabled() const { return m_enabled; }

      // @returns a string with information about this checklist
      std::string const GetInfoString() const {
        if(!m_enabled)
          return "disabled";
        switch(m_mode) {
          case kNoCuts:
            return fmt::format("for all PDGs ({})", fmt::join(m_pdg_list, ", "));
          case k1hCuts:
            return fmt::format("in [{:.3g}, {:.3g}], for all PDGs ({})", m_min, m_max, fmt::join(m_pdg_list, ", "));
          case k2hCuts:
            return fmt::format("in [{:.3g}, {:.3g}], for ({}) dihadrons", m_min, m_max, fmt::join(m_pdg_list, ", "));
        }
        throw std::runtime_error("GetInfoString failed");
      }

      /// @brief check the checklist for single particles
      /// @param evt the pythia event
      /// @param get_val a function to get the value to be checked, given a `Pythia8::Particle`; not used if `m_mode == kNoCuts`
      /// @param must_be_final if true, particle must be "final"
      /// @returns true if all particles pass the cut; if this checklist is not enabled, returns true
      bool const Check(
          Pythia8::Event const& evt,
          std::function<double(Pythia8::Particle const&)> const& get_val = [](Pythia8::Particle const& par){ return 0; },
          bool const& must_be_final = true) const
      {
        // return true if not enabled
        if(!m_enabled)
          return true;
        // initialize checklist
        if(enable_verbose_mode) fmt::println("CHECKLIST for {}", m_opt_name);
        std::vector<std::pair<int, bool>> check_list;
        for(auto const& pdg : m_pdg_list)
          check_list.emplace_back(pdg, false);
        // loop over event particles, and check the checkboxes
        std::size_t num_found = 0;
        for(auto const& par : evt) { // loop over event particles
          if(!must_be_final || (must_be_final && par.isFinal())) { // particle must be final, if `must_be_final==true`
            for(auto& [pdg, found] : check_list) { // loop over checklist
              if(!found && pdg == par.id()) { // if we haven't found this one yet, and the checklist PDG == particle PDG
                switch(m_mode) {
                  case kNoCuts:
                    {
                      // no value comparison, just check the box
                      found = true;
                      if(enable_verbose_mode) fmt::println("  [x] {} at idx {}", pdg, par.index());
                      break;
                    }
                  case k1hCuts:
                    {
                      // check the box if val is in range (m_min, m_max)
                      auto val = get_val(par);
                      found = m_min <= val && val <= m_max;
                      if(enable_verbose_mode) fmt::println("  [x] {} at idx {}, value = {}", pdg, par.index(), val);
                      break;
                    }
                  case k2hCuts:
                    {
                      throw std::runtime_error("called single-particle `Check` but mode is dihadron");
                      break;
                    }
                }
                if(found) {
                  num_found++;
                  break; // break the loop over checklist
                }
              }
            }
          }
          // if all checkboxes are checked, stop looping over event particles and return true
          if(num_found == check_list.size()) {
            return true;
          }
        }
        // not all checkboxes were checked, return false
        return false;
      }

      /// @brief check if there's a dihadron which satisfies the cut
      /// @param evt the pythia event
      /// @param dih_kin a list of `DihadronKin` objects for this event
      /// @param get_val a function to get the value to be checked, given 2 `Pythia8::Particle`s
      /// @returns true if a dihadron which satisfies the cut is found; if this checklist is not enabled, returns true
      bool const Check(
          Pythia8::Event const& evt,
          std::vector<DihadronKin> const& dih_kin,
          std::function<double(Pythia8::Particle const&, Pythia8::Particle const&)> const& get_val
          ) const
      {
        // return true if not enabled
        if(!m_enabled)
          return true;
        // must be dihadron mode
        if(m_mode != k2hCuts)
          throw std::runtime_error("called dihadron `Check` but mode is not dihadron");
        // loop over dihadrons
        if(enable_verbose_mode) fmt::println("CHECKLIST for {}", m_opt_name);
        for(auto const& dih : dih_kin) {
          // if this is a dihadron the user wants
          if( (dih.pdgA == m_pdg_list.at(0) && dih.pdgB == m_pdg_list.at(1)) ||
              (dih.pdgA == m_pdg_list.at(1) && dih.pdgB == m_pdg_list.at(0)) )
          {
            // check the cuts
            auto const& parA = evt.at(dih.idxA);
            auto const& parB = evt.at(dih.idxB);
            auto val = get_val(parA, parB);
            if(m_min <= val && val <= m_max) {
              if(enable_verbose_mode) fmt::println("  [x] ({}, {}) at idxs ({}, {}), value = {}", parA.id(), parB.id(), dih.idxA, dih.idxB, val);
              return true; // satisfactory dihadron
            }
          }
        }
        return false; // no satisfactory dihadron
      }

      /// @returns a const reference to the underlying PDG list
      std::vector<int> const& GetPDGList() const { return m_pdg_list; }

    private:
      /// the name of the checklist should match the CLI option name
      std::string const m_opt_name;
      /// the mode for this checklist
      Mode const m_mode;
      /// whether this checklist is enabled
      bool m_enabled;
      /// list of PDGs in this checklist
      std::vector<int> m_pdg_list;
      /// cut minimum
      double m_min;
      /// cut maximum
      double m_max;
  };
}
