#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'

options = {
  num_events: 10000,
  outdir:     nil,
  jobs:       1,
  slurm:      false,
  no_gemc:    false,
  sets:       []    # <-- collect all "--set PARAM=VALUE" here
}

parser = OptionParser.new do |opts|
  opts.banner = "Usage: #{$PROGRAM_NAME} [options]"

  opts.on("-n", "--num-events N", Integer,
          "Number of events (default #{options[:num_events]})") do |n|
    options[:num_events] = n
  end

  opts.on("-o", "--outdir DIR", String,
          "Output directory under /volatile/clas12/users/gmat/osg/clas-stringspinner/ (required)") do |d|
    options[:outdir] = d
  end

  opts.on("-j", "--jobs J", Integer,
          "Number of jobs to run for each spin (default #{options[:jobs]})") do |j|
    options[:jobs] = j
  end

  opts.on("--slurm", "Submit each job via SLURM instead of running locally") do
    options[:slurm] = true
  end

  opts.on("--no-gemc", "Skip GEMC, recon-util, hipo filter, and hipo2tree steps") do
    options[:no_gemc] = true
  end

  # This line allows multiple --set PARAM=VALUE flags.
  opts.on("--set PARAM=VALUE", String,
          "Pass --set PARAM=VALUE to clas-stringspinner; may be used multiple times") do |sv|
    # sv will be the literal "StringSpinner:GLGT=100", etc.
    options[:sets] << sv
  end

  opts.on("-h", "--help", "Show this help") do
    puts opts
    exit
  end
end

begin
  parser.parse!
rescue OptionParser::InvalidOption => e
  warn e.message
  puts parser
  exit 1
end

if options[:outdir].nil?
  warn "Error: --outdir is required"
  puts parser
  exit 1
end

# Prepare base directory
base_outdir = File.join('/volatile/clas12/users/gmat/osg/clas-stringspinner/', options[:outdir])
if Dir.exist?(base_outdir)
  print "Directory #{base_outdir} exists. Overwrite? [y/N]: "
  answer = STDIN.gets.chomp
  unless answer.downcase.start_with?('y')
    puts "Aborting."
    exit 1
  end
  FileUtils.rm_rf(base_outdir)
end

# Create subdirectories
%w[ err log slurm ].each do |d|
  FileUtils.mkdir_p(File.join(base_outdir, d))
end

# Locate recon-util
coatjava = ENV['COATJAVA']
abort "COATJAVA not set" if coatjava.nil? || coatjava.strip.empty?
recon_util = File.join(coatjava, 'bin', 'recon-util')
abort "recon-util not found or not executable" unless File.executable?(recon_util)

pol = 'LU'
submission_ids = []

if options[:slurm]
  # Generate and submit SLURM scripts
  ['p', 'n'].each do |spin|
    options[:jobs].times do |i|
      suffix      = format('%05d', i + 1)
      seed        = rand(1_000_000_000)
      job_name    = "#{pol}_#{spin}_#{suffix}"
      dat_file    = File.join(base_outdir, "clas_#{job_name}.txt")
      hipo_file   = File.join(base_outdir, "clas_#{job_name}.hipo")
      cooked_tmp  = File.join(base_outdir, "clas_#{job_name}_cooked_.hipo")
      cooked_file = File.join(base_outdir, "clas_#{job_name}_cooked.hipo")
      hipo2tree   = File.join(base_outdir, "clas_#{job_name}_cooked.root")
      lund2tree   = File.join(base_outdir, "clas_#{job_name}_lund.root")
      script_path = File.join(base_outdir, 'slurm', "job_#{job_name}.slurm")

      File.open(script_path, 'w') do |f|
        f.puts "#!/bin/bash"
        f.puts "#SBATCH --account=clas12"
        f.puts "#SBATCH --partition=production"
        f.puts "#SBATCH --mem-per-cpu=2000"
        f.puts "#SBATCH --job-name=#{job_name}"
        f.puts "#SBATCH --cpus-per-task=2"
        f.puts "#SBATCH --time=24:00:00"
        f.puts "#SBATCH --output=#{File.join(base_outdir, 'log', "#{job_name}.log")}"
        f.puts "#SBATCH --error=#{File.join(base_outdir, 'err', "#{job_name}.err")}"
        f.puts
        # 1) generate LUND dat
        f.puts "build/clas-stringspinner \\"
        f.puts "  --num-events #{options[:num_events]} \\"
        f.puts "  --pol-type #{pol} \\"
        f.puts "  --beam-spin #{spin} \\"
        f.puts "  --cut-family-inclusive 211,-211,111 \\"
        f.puts "  --cut-inclusive 11 \\"
        f.puts "  --cut-theta 5,35,11 \\"
        f.puts "  --cut-family-theta 5,35,211,-211,22,22 \\"
        f.puts "  --seed #{seed} \\"
        # Insert any --set flags here (each ends with a backslash)
        options[:sets].each do |sv|
          # sv is something like "StringSpinner:GLGT=100"
          f.puts "  --set #{sv} \\"
        end
        # Finally, write out-file without backslash
        f.puts "  --out-file #{dat_file}"
        f.puts

        unless options[:no_gemc]
          # 2) run GEMC
          f.puts "gemc etc/rga_fall2018.gcard \\"
          f.puts "  -SAVE_ALL_MOTHERS=1 -SKIPREJECTEDHITS=1 \\"
          f.puts "  -NGENP=#{options[:num_events]} \\"
          f.puts "  -INTEGRATEDRAW=\"*\" \\"
          f.puts "  -USE_GUI=0 \\"
          f.puts "  -RUNNO=11 \\"
          f.puts "  -INPUT_GEN_FILE=\"LUND, #{dat_file}\" \\"
          f.puts "  -OUTPUT=\"hipo, #{hipo_file}\""
          f.puts
          # 3) reconstruction
          f.puts "#{recon_util} -i #{hipo_file} -o #{cooked_tmp} -y etc/rga_fall2018.yaml"
          f.puts
          # 4) filter banks
          f.puts "hipo-utils -filter -b 'REC::Calorimeter,REC::Particle,MC::Particle,RUN::config' #{cooked_tmp} -o #{cooked_file}"
          f.puts
          # 5) convert hipo to ROOT
          f.puts "./build/hipo2tree #{cooked_file} #{hipo2tree}"
          f.puts
        end

        # 6) always convert LUND to ROOT
        f.puts "./build/lund2tree #{dat_file} #{lund2tree}"
        f.puts
        # 7) clean up intermediate files
        if !options[:no_gemc]
          f.puts "rm #{cooked_tmp} #{hipo_file} #{dat_file} #{cooked_file}"
          f.puts
        end
        # 8) completion message
        f.puts "echo \"Completed job #{job_name}\""
      end

      print "Submitting #{script_path}..."
      result = `sbatch #{script_path}`
      if result =~ /Submitted batch job (\d+)/
        submission_ids << $1
        puts " #{result.strip}"
      else
        abort("sbatch failed for job #{job_name}")
      end
    end
  end

  # Submit merge job after all
  if submission_ids.any?
    dep = submission_ids.join(':')
    merge_script = File.join(base_outdir, 'slurm', "merge_all.slurm")
    File.open(merge_script, 'w') do |f|
      f.puts "#!/bin/bash"
      f.puts "#SBATCH --account=clas12"
      f.puts "#SBATCH --partition=production"
      f.puts "#SBATCH --dependency=afterok:#{dep}"
      f.puts "#SBATCH --mem-per-cpu=4000"
      f.puts "#SBATCH --cpus-per-task=4"
      f.puts "#SBATCH --output=#{File.join(base_outdir, 'log', 'merge_all.log')}"
      f.puts "#SBATCH --error=#{File.join(base_outdir, 'err', 'merge_all.err')}"
      f.puts
      f.puts "hadd #{File.join(base_outdir, "lund_merged.root")} #{File.join(base_outdir, 'clas_*_lund.root')}"
      unless options[:no_gemc]
        f.puts "hadd #{File.join(base_outdir, "cooked_merged.root")} #{File.join(base_outdir, 'clas_*_cooked.root')}"
      end
      f.puts "echo \"Merged ROOT files into #{options[:outdir]}_*_merged.root\""
    end
    puts "Submitting merge job..."
    system("sbatch #{merge_script}") or abort("sbatch failed for merge job")
  end

else
  # Run locally
  ['p', 'n'].each do |spin|
    options[:jobs].times do |i|
      suffix      = format('%05d', i + 1)
      seed        = rand(1_000_000_000)
      job_name    = "#{pol}_#{spin}_#{suffix}"
      dat_file    = File.join(base_outdir, "clas_#{job_name}.txt")
      hipo_file   = File.join(base_outdir, "clas_#{job_name}.hipo")
      cooked_tmp  = File.join(base_outdir, "clas_#{job_name}_cooked_.hipo")
      cooked_file = File.join(base_outdir, "clas_#{job_name}_cooked.hipo")
      hipo2tree   = File.join(base_outdir, "clas_#{job_name}_cooked.root")
      lund2tree   = File.join(base_outdir, "clas_#{job_name}_lund.root")

      puts "Job #{i+1}/#{options[:jobs]} [spin=#{spin}] (seed=#{seed})"

      # 1) generate LUND dat with any --set flags appended
      cmd_parts = [
        "build/clas-stringspinner",
        "--num-events #{options[:num_events]}",
        "--pol-type #{pol}",
        "--beam-spin #{spin}",
        "--cut-inclusive 11",
        "--cut-family-inclusive 211,-211,111",
        "--cut-theta 5,35,11",
        "--cut-family-theta 5,35,211,-211,22,22",
        "--seed #{seed}",
        "--out-file #{dat_file}"
      ]

      # Add each "--set PARAM=VALUE" exactly as passed
      options[:sets].each do |sv|
        cmd_parts << "--set #{sv}"
      end

      cmd1 = cmd_parts.join(' ')
      system(cmd1) or abort("clas-stringspinner failed on job #{job_name}")

      unless options[:no_gemc]
        # 2) run GEMC
        cmd2 = [
          "gemc etc/rga_fall2018.gcard",
          "-SAVE_ALL_MOTHERS=1",
          "-SKIPREJECTEDHITS=1",
          "-NGENP=#{options[:num_events]}",
          "-INTEGRATEDRAW=\"*\"",
          "-USE_GUI=0",
          "-RUNNO=11",
          "-INPUT_GEN_FILE=\"LUND, #{dat_file}\"",
          "-OUTPUT=\"hipo, #{hipo_file}\""
        ].join(' ')
        system(cmd2) or abort("gemc failed on job #{job_name}")

        # 3) reconstruction
        system("#{recon_util} -i #{hipo_file} -o #{cooked_tmp} -y etc/rga_fall2018.yaml") or abort("recon-util failed on job #{job_name}")

        # 4) filter banks
        system("hipo-utils -filter -b 'REC::Calorimeter,REC::Particle,MC::Particle,RUN::config' #{cooked_tmp} -o #{cooked_file}") or abort("hipo-utils failed on job #{job_name}")

        # 5) convert hipo to ROOT
        system("./build/hipo2tree #{cooked_file} #{hipo2tree}") or abort("hipo2tree failed on job #{job_name}")
      end

      # 6) always convert LUND to ROOT
      system("./build/lund2tree #{dat_file} #{lund2tree}") or abort("lund2tree failed on job #{job_name}")

      puts "Completed job #{job_name}"
    end
  end

  # Local merge
  puts "Merging LUND ROOT files..."
  merged_lund = File.join(base_outdir, "lund_merged.root")
  system("hadd #{merged_lund} #{File.join(base_outdir, 'clas_*_lund.root')}") or abort("hadd failed for LUND files")
  unless options[:no_gemc]
    puts "Merging cooked ROOT files..."
    merged_cooked = File.join(base_outdir, "cooked_merged.root")
    system("hadd #{merged_cooked} #{File.join(base_outdir, 'clas_*_cooked.root')}") or abort("hadd failed for cooked files")
  end
  puts "Merge complete."
end
