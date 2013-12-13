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
		options[:path_to_alignment] = nil
		options[:path_to_genestruct] = nil
		options[:path_to_output] = "genepainter"
		options[:path_to_log] = "genepainter.log"
		options[:output_format] = ["txt_simple"]

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

# the old options
# -pdb file.pdb [chain]         Two scripts for execution in PyMol will be provided. In color_exons.py the
#                               consensus exons are colored and in color_splicesites.py the splice junctions
#                               of the consensus exons are marked for <chain> (default: chain A)
# -pdb_prot prot_name           Use protein <prot_name> as reference for alignment with the pdb sequence
#                               Default: First protein in <alignment>
# -f                            Force alignment between pdb and first protein sequence of the MSA or protein <prot_name> (if specified)
#                               This ignores the default that intron positions will only be mapped if the alignment score > 70%
# -consensus value              Color only intron positions conserved in <value> percent of all genes (default: 80%)
# -ref_prot_struct              Color only the intron positions occurring in the gene of the reference protein.
#                               May be combined with "-consensus"
# -penalize_endgaps             Penalize gaps at the end of the alignment (behaves like the standard Needleman-Wunsch algorithm)
#                               Default: gaps at the end of the alignment are not penalized'


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
				options[:output_format] << "txt_intron_phases"
			end
			opts.on("--phylo", 
				"Mark exons by '0' and introns by '1'") do
				options[:output_format] << "txt_phylo"
			end
			opts.on("--spaces",
			 "Mark exons by space (' ') instead of '-'" ) do
				options[:output_format] << "txt_only_introns"
			end
			opts.on("--no-standard-output", 
				"Specify to skip standard output format." ) do
				options[:output_format].delete("txt_simple")
			end

			opts.on("--alignment",
				"Output the alignment file with additional lines containing intron phases") do
				options[:output_format] << "alignment_with_intron_phases"
			end

			opts.separator ""
			opts.separator "Graphical output format:"
			opts.on("--svg H,W", Array,
				"Drawn SVG of size height x width") do |list|
				list = list.map(&:to_i)
				if list.size != 2 || list.inject(:*) == 0 then
					# number of args wrong or at least one is zero
					Helper.abort "Invalid argument: --svg expects two numbers"
				end
				options[:output_format] << "svg"
				if options[:svg_options] then
					options[:svg_options][:size][:height] = list[0]
					options[:svg_options][:size][:width] = list[1]
				else
					options[:svg_options] = { size: {height: list[0], width: list[1]} }
				end
			end
			svg_formats = ["normal", "reduced"]
			opts.on("--svg-format FORMAT", svg_formats,
				"	FORMAT: #{svg_formats}",
				"	'normal' draws details of aligned exons and introns [default]",
				"	'reduced' focuses on common introns only") do |f|
				if f == "reduced" then
					if options[:svg_options] then
						options[:svg_options][:reduced] = true
					else
						options[:svg_options] = { reduced: true }
					end
				else
					# nothing to do, as 'normal' is default format and already listed as output format
				end
			end

			opts.separator ""
			opts.separator "Meta information and statistics:"
			opts.on("--consensus N", Float, 
				"Introns conserved in N % genes.", 
				"	Specify N as decimal number between 0 and 1") do |n|
				options[:consensus] = n
			end
			opts.on("--merge", "Merge all introns into a single exon intron pattern") do 
				options[:merge] = true 
			end

			opts.separator ""
			opts.separator "General options:"
			opts.on("-o", "--outfile <file_name>", 
				"Name of the output file") do |file|
				options[:path_to_output] = File.basename(file, ".*")
				options[:path_to_log] = options[:path_to_output] + ".log"
			end

			opts.on("--range <start-stop[,start-stop,...]>", Array,
				"Comma-separated list (without blanks) of alignment ranges.",
				"	Keyword 'end' for last alignment position can be used") do |range_list|

				parsed_list = []
				last_r_stop = nil
				# sorting the range_list according to starting positions as kind of user-service
				range_list.sort.each do |range|
					range = range.strip
					splitted_parts = range.split("-")
					if splitted_parts.size != 2 then
						# invalid syntax used for range definition: cannot parse it
						Helper.abort("Invalid range definition in #{range}: start and stop must be separated by '-'")
					else
						# syntax was ok, check for content
						r_start = splitted_parts[0].to_i - 1 # convert from human to ruby counting
						# set r_stop to Infinity if end of range should cover alignment till end
						if range.match(/end/) then
							r_stop = Float::INFINITY
						else
							r_stop = splitted_parts[1].to_i - 1 # convert from human to ruby counting
						end

						if r_start >= r_stop && r_stop != -1 then
							# invalid syntax, start must be before stop
							Helper.abort("Invalid range definition in #{range}: start must be lower than stop")
						elsif last_r_stop && r_start <= last_r_stop 
							# invalid syntax: ranges are overlapping
							Helper.abort("Invalid range definition in #{range_list}: overlapping ranges")
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
		elsif options[:output_format].empty?
			Helper.abort "Select at least one output format"
				
		end

		return options

		# use the own format of fatal error messages!				
		rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption => exc
			Helper.abort exc
	end # parse()

end
