require 'optparse'

class OptParser
	def self.parse(args)

		options = Hash.new

		# initialize default params
		options[:path_to_alignment] = nil
		options[:path_to_genestruct] = nil
		options[:path_to_output] = "genepainter"
		options[:path_to_log] = "genepainter.log"
		options[:output_format] = "txt_simple"

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
			opts.separator ""

			opts.on("-i <path_to_alignment>", 
				"Path to fasta-formatted multiple sequence alignment") do |file|
				options[:path_to_alignment] = file
				Helper.file_exist_or_die(file)
			end

			opts.on("-p <path_to_genestructures>",
				"Path to folder containing gene structures in YAML or GFF format") do |path|
				options[:path_to_genestruct] = path
				Helper.dir_exist_or_die(path)
			end

			opts.separator ""
			opts.separator "Options:"
			opts.on("-o <file_name>", 
				"Name of the output file") do |file|
				options[:path_to_output] = File.basename(file, ".*")
				options[:path_to_log] = options[:path_to_output] + ".log"
			end
			opts.on("-t [normal|phase|space|phylo]", [:normal, :phase, :space, :phylo],
				"Format of text output, default: '-t normal'",
				"	<normal>	Mark exons by '-' and introns by '|'", 
				"	<phase> 	Mark introns by their phase instead of '|'",
				"	<space>		Mark exons by ' ' instead of '-'",
				"	<phylo>		Mark exons by '0' and introns by '1'") do |opt|

				case opt 
				when :normal
					str = "txt_simple"
				when :phase
					str = "txt_intron_phases"
				when :space
					str = "txt_only_introns"
				when :phylo
					str = "txt_phylo"
				end
					
				options[:output_format] = str 
			end
			opts.on("-a", 
				"Output the alignment file with additional lines containing intron phases") do |a|
				options[:output_format] = "alignment_with_intron_phases"
			end
			#TODO
			# opts.on("-svg")
# Options:
	# when "svg"
	# 	svg_width = args.shift.to_i 
	# 	svg_height = args.shift.to_i
	# 	# determine output format
	# 	case args.shift
	# 	when "normal"
	# 		output_format = "svg_focus_on_exons_and_introns"
	# 	when "reduced"
	# 		output_format = "svg_focus_on_common_introns"
	# 	else
	# 		puts "Unknown specification for SVG output. Using \"normal\" instead."
	# 		output_format = "svg_focus_on_exons_and_introns"
	# 	end	
	# 	# check if SVG size is valid
	# 	if svg_width == 0 || svg_height == 0 then
	# 		puts "Illegal specification of SVG size. Creating SVG of size 1000 x 1000 instead."
	# 		svg_width = 1000
	# 		svg_height = 1000
	# 	end
	# end

# -svg <width> <height> [normal|reduced] \tCreate and SVG file of size <width> x <height>
# \t\t\t\tDraw exons and introns with detail with option <normal>; Focus on common introns with option <reduced>


			opts.separator ""
			opts.on_tail("-h", "--help", "Show this message") do 
				puts opts
				exit
			end

		end # optionparser

		opt_parser.parse(args)

		# make sure that mandatory arguments are present
		if ! options[:path_to_alignment] then
			Helper.abort "Mandatory argument '-i' is missing"
		elsif ! options[:path_to_genestruct]
			Helper.abort "Mandatory argument '-p' is missing"
		end

		return options

		# use the own format of fatal error messages!				
		rescue OptionParser::MissingArgument => exc
			Helper.abort exc
		rescue OptionParser::InvalidOption => exc 
			Helper.abort exc
	end # parse()

end
