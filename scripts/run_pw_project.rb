#!/usr/bin/env ruby
require 'json'

# map each pion_pair => its run_version & beam polarization
projects = {
  "piplus_pi0"  => { run_version: "Fall2018Spring2019_RGA_inbending", beam_polarization: 0.8516 },
  "piplus_piplus"  => { run_version: "Fall2018Spring2019_RGA_inbending", beam_polarization: 0.8516 },
  "piplus_piminus"  => { run_version: "Fall2018Spring2019_RGA_inbending", beam_polarization: 0.8516 },
  "piminus_pi0" => { run_version: "Fall2018_RGA_outbending",                 beam_polarization: 0.8922 },
  "piminus_piminus" => { run_version: "Fall2018_RGA_outbending",                 beam_polarization: 0.8922 }
}



# your three 10-bin definitions
binning_list = [
  { "x"  => [0.06,0.11,0.13,0.15,0.17,0.2,0.22,0.26,0.33,0.5,0.75] },
  { "z"  => [0.2,0.34,0.38,0.42,0.46,0.5,0.55,0.6,0.68,0.75,0.95] },
  { "Mh" => [0.16,0.39,0.47,0.54,0.6,0.68,0.75,0.84,0.97,1.4,2.53] }
]


date_str = Time.now.strftime("%m_%d_%Y")
subdir   = "asymmetry___pw_acc_is_rho"
generator = File.expand_path("generate_slurm.rb", __dir__)
file_name = "/volatile/clas12/users/gmat/osg/clas-stringspinner/zero_quark_masses_10M/lund_merged.root"
global_cuts = { "cut1" => ":is_rho:==1" }   
global_cuts_json = global_cuts.to_json



projects.each do |pion_pair, cfg|
  binning_list.each do |bin|
    json_bin = bin.to_json

    tree_name = "dih_#{pion_pair}_acceptance"
    cmd = [
      "ruby", generator,
      "--date",          date_str,
      "--subdir",        subdir,
      "--pion_pair",     pion_pair,
      "--run_version",   cfg[:run_version],
      "--beam_polarization", cfg[:beam_polarization].to_s,
      "--binning_info",  "'#{json_bin}'",
      "--tree_name", tree_name,
      "--file_name", file_name,
      "--global_cuts", "'#{global_cuts_json}'",
      "--no-slurm",

    ].join(" ")

    puts "\n>> Launching for #{pion_pair}, bin #{bin.keys.first}:"
    puts cmd
    system(cmd)
  end
end
