#!/usr/bin/env ruby
require 'fileutils'
require 'json'
require 'optparse'

options = {
  file_name:         '',
  no_slurm:          false,
  overwrite:         false,
  date:              nil,
  subdir:            nil,
  pion_pair:         nil,
  run_version:       nil,
  beam_polarization: nil,
  binning_info:      nil,
  tree_name:         nil,
  global_cuts: nil,

}

OptionParser.new do |opts|
  opts.banner = "Usage: generate_slurm.rb [options]"

  opts.on("--file_name=FILE", "-f",               "Overwrite file name") do |v|
    options[:file_name] = v
  end
    
  opts.on("--no-slurm",               "Run jobs directly instead of generating/submitting slurm scripts") do
    options[:no_slurm] = true
  end

  opts.on("--overwrite", "-o",        "Automatically overwrite existing directories") do
    options[:overwrite] = true
  end

  opts.on("--date=MM_DD_YYYY", "-d",  "Date string (default: today)") do |v|
    options[:date] = v
  end

  opts.on("--subdir=NAME", "-s",      "Additional subdirectory under the date") do |v|
    options[:subdir] = v
  end

  opts.on("--pion_pair=PAIR",         "Pion pair (overrides default)") do |v|
    options[:pion_pair] = v
  end

  opts.on("--run_version=VER",        "Run version / dataset (overrides default)") do |v|
    options[:run_version] = v
  end

  opts.on("--beam_polarization=POL", Float, "Beam polarization (overrides default)") do |v|
    options[:beam_polarization] = v
  end

  opts.on("--binning_info=JSON",      "Binning info as JSON, e.g. '{\"x\":[0.06,0.11,…]}'") do |v|
    options[:binning_info] = JSON.parse(v)
  end

  opts.on("--tree_name=TREE_NAME",      "TTree name") do |v|
    options[:tree_name] = v
  end
  opts.on("--global_cuts=JSON", "Global cuts as JSON, e.g. '{\"z2_gt_z1\":\":z2:>:z1:\"}'") do |v|
    options[:global_cuts] = v
  end

  opts.on("-h","--help",              "Show this help") do
    puts opts
    exit
  end
end.parse!

date_str = options[:date] || Time.now.strftime("%m_%d_%Y")

# wrap user‐passed binning_info or fall back to default
binning_infos = if options[:binning_info]
                  [ options[:binning_info] ]
                else
                  [ { "x" => [0.06,0.11] } ]
                end

pion_pair         = options[:pion_pair]      || "piplus_pi0"
run_version       = options[:run_version]    || "Fall2018_RGA_inbending"
beam_polarization = options[:beam_polarization] || 0.8516
tree_name = options[:tree_name] || "dihadron_cuts"
binning_infos.each do |info|
  info.each do |var, edges|
    # build: out_slurm/<date>/<subdir?>/<var>/<pion_pair>
    path_components = ["out_slurm", date_str]
    path_components << options[:subdir] if options[:subdir]
    path_components += [var, pion_pair]
    base_dir  = File.join(*path_components)
    slurm_dir = File.join(base_dir, "slurm")
    yaml_dir  = File.join(base_dir, "yaml")

    # handle overwrite/prompt
    if Dir.exist?(base_dir)
      if options[:overwrite]
        puts "Overwriting existing directory #{base_dir}..."
        FileUtils.rm_rf(base_dir)
      else
        print "Directory #{base_dir} already exists. Overwrite? [y/N]: "
        if STDIN.gets.chomp.downcase == 'y'
          FileUtils.rm_rf(base_dir)
        else
          puts "Skipping #{var} (directory exists)."
          next
        end
      end
    end

    FileUtils.mkdir_p([slurm_dir, yaml_dir])

    # collect job IDs
    job_ids = []

    # ========== worker jobs ==========
    edges.each_cons(2) do |lo, hi|
      bin_label = "#{var}_#{lo}-#{hi}"
      json_bin  = { var => [lo, hi] }.to_json
      # ---- build the python command ------------------------------------------------
      cmd_parts = [
          "python3 /work/clas12/users/gmat/clas12/clas12_dihadrons/plotting/release_pw_pi0copy/macros/run_asymmetry.py",
          "--binning_info '#{json_bin}'",
          "--pion_pair '#{pion_pair}'",
          "--run_version '#{run_version}'",
          "--beam_polarization #{beam_polarization}",
          "--output_dir '#{yaml_dir}'",
          "--subtitle '#{bin_label}'",
          "--tree_name '#{tree_name}'",
          "--file_path #{options[:file_name]}",
          "--ignore-sideband True"
        ]
        
      # add the override only when the user supplied it
      cmd_parts << "--global_cuts '#{options[:global_cuts]}'" if options[:global_cuts]
        
      cmd = cmd_parts.join(' ')

      if options[:no_slurm]
        puts "=== running worker for #{bin_label} ==="
        system(cmd)
      else
        slurm_file = File.join(slurm_dir, "#{bin_label}.slurm")
        File.open(slurm_file, 'w', 0755) do |f|
          f.puts "#!/bin/bash"
          f.puts "#SBATCH --job-name=#{bin_label}"
          f.puts "#SBATCH --output=#{slurm_dir}/#{bin_label}.out"
          f.puts "#SBATCH --error=#{slurm_dir}/#{bin_label}.err"
          f.puts "#SBATCH --partition=production"
          f.puts "#SBATCH --account=clas12"
          f.puts "#SBATCH --cpus-per-task=4"
          f.puts "#SBATCH --mem-per-cpu=4000"
          f.puts "#SBATCH --time=24:00:00"
          f.puts
          f.puts cmd
          f.puts
        end

        # submit and capture job ID
        sbatch_out = `sbatch #{slurm_file}`.strip
        if sbatch_out =~ /Submitted batch job (\d+)/
          job_id = $1
          job_ids << job_id
          puts "Submitted worker #{bin_label} as job #{job_id}"
        else
          warn "Failed to submit #{slurm_file}: #{sbatch_out}"
        end
      end
    end

    # ================== collector job ==================
    collector_sh = File.join(slurm_dir, "collector_#{var}.slurm")
    collector_bash = <<~BASH
      #!/bin/bash
      #SBATCH --job-name=collect_#{pion_pair}_#{run_version}_#{var}
      #SBATCH --output=#{slurm_dir}/collect_#{var}.out
      #SBATCH --error=#{slurm_dir}/collect_#{var}.err
      #SBATCH --partition=production
      #SBATCH --account=clas12
      #SBATCH --cpus-per-task=1
      #SBATCH --time=01:00:00

      python3 /work/clas12/users/gmat/clas12/clas12_dihadrons/plotting/release_pw_pi0copy/macros/merge_asymmetry_yaml.py \\
        --yaml_base "#{yaml_dir}" \\
        --variable "#{var}" \\
        --pion_pair "#{pion_pair}" \\
        --run_version "#{run_version}"
    BASH

    File.write(collector_sh, collector_bash)
    FileUtils.chmod("+x", collector_sh)

    if options[:no_slurm]
      puts "=== running collector for #{var} locally ==="
      system("bash #{collector_sh}")
    else
      if job_ids.any?
        dep = job_ids.join(':')
        submit_cmd = "sbatch --dependency=afterok:#{dep} #{collector_sh}"
        dep_out = `#{submit_cmd}`.strip
        if dep_out =~ /Submitted batch job (\d+)/
          puts "Submitted collector for #{var} (depends on #{dep}) as job #{$1}"
        else
          warn "Failed to submit collector: #{dep_out}"
        end
      else
        warn "No worker jobs to depend on; submitting collector immediately."
        system("sbatch #{collector_sh}")
      end
    end
  end
end
