require 'optparse'

class OptParser

	def self.backwards_compability_check(args)
		old_options = ["-a", "-n", "-phylo", "-s", "-svg", "-start", "-stop", 
				"-pdb", "-pdb_prot", "-f", "-consensus", "-ref_prot_struct", "-penalize_endgaps"
		]
		# overlap between old options and args specified 
		specified_old_options = old_options & args
		if specified_old_options.any? then
			# print an error message
			Helper.warn "Parameter(s) #{specified_old_options.join(", ")} have changed. Plese refer to their new definition."
			Helper.warn "---"
			parse(["-h"]) # enforce that (new) options definitions are printed
		end
	end

	def self.parse(args)

		options = Hash.new

		# initialize default params

		# mandatory parameters
		options[:path_to_alignment] = nil # input: alignment
		options[:path_to_genestruct] = nil # input: gene structures

		# optional parameters
		options[:path_to_output] = "genepainter" # output file
		options[:path_to_log] = "genepainter.log" # log file

		options[:output_format_list] = ["txt_simple"] # which output should be generated
		options[:range] = nil # restrict input alignment: use only "columns" within range
		options[:ignore_common_gaps] = false # restrict input alignment: ignore gaps common to all sequences
		options[:consensus] = false # add consenus profile to output ?
		options[:merge] = false # add merged profile to output ?

		options[:svg_options] = {} # options for svg output, only set if svg_output is requested
		options[:pdb_options] = {} # options for pdb output, only set if pdb_output is requested

		# a list of all possible output formats
		# possible_output_formats = [
		# 	"txt_simple", # plain text, exons => "-", introns => "|"
		# 	"txt_intron_phases", # plain text, exons => "-", introns => "[0|1|2|?]"
		# 	"txt_only_introns", # plain text, exons => " ", introns => "|"
		# 	"txt_phylo", # plain text, exons => "0", introns => "1"
		# 	"alignment_with_intron_phases", # fasta-formatted alignment with additional for every sequence containing intron phases
		# 	"svg_focus_on_exons_and_introns", # svg, aligned exons and introns shown
		# 	"svg_focus_on_common_introns" # svg, conserved introns drawn in same colour
		# ]

		opt_parser = OptionParser.new do |opts|

			opts.banner = "\nGenePainter v.1.2 maps gene structures onto multiple sequence alignments"
			opts.separator "Contact: Martin Kollmar (mako[at]nmr.mpibpc.mpg.de)"

			opts.separator ""
			opts.separator "Usage: gene_painter.rb -i path_to_alignment -p path_to_genestructure_folder [options]"
			opts.separator "Standard output format: Mark exons by '-' and introns by '|'"
			opts.separator ""

			opts.on("-i", "--input <path_to_alignment>",
				"Path to fasta-formatted multiple sequence alignment") do |file|
				options[:path_to_alignment] = file
				Helper.file_exist_or_die(file)
			end

			opts.on("-p", "--path <path_to_genestructures>",
				"Path to folder containing gene structures in YAML or GFF format") do |path|
				options[:path_to_genestruct] = path
				Helper.dir_exist_or_die(path)
			end

			opts.separator ""
			opts.separator "Options:"

			opts.separator "Text-based output format:"
			opts.on("--intron-phase", 
				"Mark introns by their phase instead of '|'" ) do
				options[:output_format_list] << "txt_intron_phases"
			end
			opts.on("--phylo", 
				"Mark exons by '0' and introns by '1'") do
				options[:output_format_list] << "txt_phylo"
			end
			opts.on("--spaces",
			 "Mark exons by space (' ') instead of '-'" ) do
				options[:output_format_list] << "txt_only_introns"
			end
			opts.on("--no-standard-output", 
				"Specify to skip standard output format." ) do
				options[:output_format_list].delete("txt_simple")
			end

			opts.on("--alignment",
				"Output the alignment file with additional lines containing intron phases") do
				options[:output_format_list] << "alignment_with_intron_phases"
			end

			opts.separator ""
			opts.separator "Graphical output format:"
			opts.on("--svg H,W", Array,
				"Drawn SVG of size height x width") do |list|
				list = list.map(&:to_i)
				
				if list.size != 2 || list.inject(:*) == 0 then
					# number of args wrong or at least one is zero
					Helper.abort "Invalid argument: --svg expects two comma-separated numbers without spaces"
				end
				options[:output_format_list] << "svg"

				vivify_hash(options, :svg_options, :size, {height: list[0], width: list[1]} )
			end
			svg_formats = ["normal", "reduced"]
			opts.on("--svg-format FORMAT", svg_formats,
				"FORMAT: #{svg_formats}",
				"'normal' draws details of aligned exons and introns [default]",
				"'reduced' focuses on common introns only") do |f|
				if f == "reduced" then
					vivify_hash(options, :svg_options, :reduced, true)
				else
					# nothing to do, as 'normal' is default format and already listed as output format
				end
			end

			opts.separator ""
			opts.on("--pdb FILE", String,
				"Mark consensus or merged gene structure in pdb FILE", 
				"Consenus gene structure contains introns conserved in N % of all genes",
				"Specify N with option --consensus N; [default: 80%]",
				"Two scripts for execution in PyMol are provided:", 
				"'color_exons.py' to mark consensus exons",
				"'color_splicesites.py' to mark splice junctions of consensus exons") do |file|
				options[:output_format_list] << "pdb"
				vivify_hash(options, :pdb_options, :path_to_pdb, file)
				Helper.file_exist_or_die(file)
			end
			opts.on("--pdb-chain CHAIN", String,
				"Mark gene structures for chain CHAIN",
				"[Default: Use chain A]") do |opt|
				vivify_hash(options, :pdb_options, :pdb_chain, opt)
			end
			opts.on("--pdb-ref-prot PROT", String,
				"Use protein PROT as reference for alignment with pdb sequence", 
				"[Default: First protein in alignment]") do |n|
				vivify_hash(options, :pdb_options, :pdb_reference_protein, n)
			end
# opts.on("--pdb-force-alignment",
#        "Force alignment between pdb and reference protein", 
#        "Without this options set, alignment score > 70% is required for mapping." ) do |opt|
# 	vivify_hash(options, :pdb, :pdb_force_alignment, true)
# end
			opts.on("--pdb-ref-prot-struct",
				"Color only intron positions occuring in the reference protein structure.") do |opt|
				vivify_hash(options, :pdb_options, :pdb_ref_prot_struct_only, true)
			end
# opts.on("--pdb-penalize-endgaps", 
# 	"Penalize gaps at end of the alignment (standard Needleman-Wunsch algorithm)", 
# 	"[Default: Gaps at alignment ends are not penalized.]") do |opt|
# 	vivify_hash(options, :pdb, :pdb_penalize_endgaps, true)
# end

			opts.separator ""
			opts.separator "Meta information and statistics:"
			opts.on("--consensus N", Float, 
				"Introns conserved in N % genes.", 
				"Specify N as decimal number between 0 and 1") do |n|
				if 0 < n && n <= 1.0 then
					options[:consensus] = n
				else
					Helper.abort("Invalid argument: --consensus expects a decimal number between 0 and 1")
				end
			end
			opts.on("--merge", "Merge all introns into a single exon intron pattern") do 
				options[:merge] = true 
			end

			opts.separator ""
			opts.separator "General options:"
			opts.on("-o", "--outfile <file_name>", String, 
				"Name of the output file") do |file|
				options[:path_to_output] = File.basename(file, ".*")
				options[:path_to_log] = options[:path_to_output] + ".log"
			end

			opts.on("--range <start-stop[,start-stop,...]>", Array,
				"Comma-separated list (without blanks) of alignment ranges of format start-stop.",
				"Keyword 'end' for last alignment position can be used") do |range_list|

				parsed_list = []
				last_r_stop = nil
				# sorting the range_list according to starting positions as kind of user-service
				range_list.sort.each do |range|
					range = range.strip
					splitted_parts = range.split("-")
					if splitted_parts.size != 2 then
						# invalid syntax used for range definition: cannot parse it
						Helper.abort("invalid range definition in #{range}: start and stop must be separated by '-'")
					else
						# syntax was ok, check for content
						r_start = splitted_parts[0].to_i - 1 # convert from human to ruby counting
						# set r_stop to Infinity if end of range should cover alignment till end
						if range.match(/end/) then
							r_stop = Float::INFINITY
						else
							r_stop = splitted_parts[1].to_i - 1 # convert from human to ruby counting
						end
						if r_start < 0 || r_stop < 0 then 
							Helper.abort("invalid range definition in #{range}: start and stop must be natural numbers")
						end
						if r_start >= r_stop then
							# invalid syntax, start must be before stop
							Helper.abort("invalid range definition in #{range}: start must be lower than stop")
						elsif last_r_stop && r_start <= last_r_stop 
							# invalid syntax: ranges are overlapping
							Helper.abort("invalid range definition in #{range_list}: overlapping ranges")
						else
							# valid syntax
							last_r_stop = r_stop
							# but: don't store r_stop as Infinity
							# this value is really cool for validity checking and really bad for index-based string access
							r_stop = -1 if r_stop == Float::INFINITY
							parsed_list << [r_start, r_stop]	
						end
					end
				end
				options[:range] = parsed_list
			end

			opts.on("--ignore-common-gaps", 
				"Ignore common gaps in alignment") do
				options[:ignore_common_gaps] = true
			end
	
			opts.separator ""
			opts.on_tail("-h", "--help", "Show this message") do 
				puts opts
				exit
			end

		end # optionparser

		## main part of function ##
		if args.empty? then
			# display help and exit if program is called without any argument
			puts opt_parser.help
			exit
		end

		opt_parser.parse(args)

		# make sure that mandatory arguments are present
		if ! options[:path_to_alignment] then
			Helper.abort "Mandatory argument '-i' is missing"
		elsif ! options[:path_to_genestruct]
			Helper.abort "Mandatory argument '-p' is missing"
		elsif options[:output_format_list].empty?
			Helper.abort "Select at least one output format"
		end

		# check dependecies
		# if pdb output should be generated, one of the following options must be set
		if options[:output_format_list].include?("pdb") then 
			if ! (
				options[:merge] || 
				options[:consensus] || 
				options[:pdb_options][:pdb_ref_prot_struct_only] ) then

			Helper.abort "Mandatory argument for --pdb is missing: Specify --consensus N, --merge or --pdb-ref-prot-struct"
			end
		end

		return options

		# use the own format of fatal error messages!				
		rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption => exc
			Helper.abort exc
	end # parse()

	def self.vivify_hash(hash, outer_key, inner_key, inner_val)
		if ! hash[outer_key] then 
			hash[outer_key] = {}
		end
		hash[outer_key][inner_key] = inner_val
	end

end
