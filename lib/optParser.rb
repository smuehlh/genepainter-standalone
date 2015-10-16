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
		options[:path_to_output] = "" # path to output files including file prefix
		options[:path_to_log] = "genepainter.log" # log file
		options[:output_basename], options[:output_basepath] = "genepainter", "" # default prefix: "genepainter", located in genepainter source directory

		options[:output_format_list] = ["txt_simple"] # which output should be generated
		options[:output_not_reduced] = false
		options[:range] = { reverse_position: [], is_delete_range: false } # restrict input alignment: use only "columns" within range
		options[:ignore_common_gaps] = true # restrict input alignment: ignore gaps common to all sequences
		options[:is_long_text_based_output] = true # text-based output contains one common gap between two introns
		options[:consensus] = false # add consenus profile to output ?
		options[:merge] = false # add merged profile to output ?

		options[:svg_options] = {} # options for svg output, only set if svg_output is requested
		options[:pdb_options] = {} # options for pdb output, only set if pdb_output is requested
		options[:tax_options] = {} # options for taxonomy output, only set if taxonomy is requested

		options[:best_intron_pos] = true # align intron positions which differ only by alignment gaps

		options[:selection] = { analyse_all_output_sel: nil, analyse_sel_output_sel: nil, analyse_sel_on_all_output_sel: nil, no_selection: nil }
		options[:select_by] = { regex: nil, list: nil, species: nil, no_selection: nil } # selection criteria: species, list, or regular expression

		# hidden parameters
		@hidden_params = ["--svg-merged", "--svg-nested", "--intron-numbers-per-taxon", "--not-reduced"]


		opt_parser = OptionParser.new do |opts|

			opts.banner = "\nGenePainter v.2.0 maps gene structures onto multiple sequence alignments"
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
				"Path to folder containing gene structures in YAML or GFF format.",
				"Required file extension is one of .yaml, .gff, .gff3") do |path|
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
			opts.on("--fuzzy N", Integer,
				"Introns at most N base pairs apart from each other are aligned") do |n_bp|
				if n_bp > 0 then 
					options[:output_format_list] << "txt_fuzzy"
					options[:fuzzy_window] = n_bp
				else
					Helper.warn "Ignoring argument: --fuzzy 0"
				end
			end
			opts.on("--not-reduced", 
				"Output patterns not reduced. This is hardly every usefull.") do
				options[:output_not_reduced] = true
			end

			opts.separator ""
			opts.separator "Graphical output format:"
			opts.on("--svg", 
				"Drawn a graphical representation of genes in SVG format", 
				"Per default, detailed representation will be produced", 
				"Use parameter '--svg-format' to get less details") do |opt|
				options[:output_format_list] << "svg"
			end
			svg_formats = ["normal", "reduced", "both"]
			opts.on("--svg-format FORMAT", svg_formats,
				"FORMAT: #{svg_formats}",
				"'normal' draws details of aligned exons and introns [default]",
				"'reduced' focuses on common introns only", 
				"'both' draws both formats") do |format|
				if format == "reduced" then
					vivify_hash(options, :svg_options, :reduced, true)
				elsif format == "both"
					vivify_hash(options, :svg_options, :both, true)
				else
					# nothing to do, as 'normal' is default format and already listed as output format
				end
			end
			opts.on("--svg-merged") do |opt|
				vivify_hash(options, :svg_options, :generate_merged_pattern, true)
			end
			opts.on("--svg-nested") do |opt|
				vivify_hash(options, :svg_options, :generate_nested_svg, true)
			end

			opts.separator ""
			opts.on("--pdb FILE", String,
				"Mark consensus or merged gene structure in pdb FILE", 
				"Consenus gene structure contains introns conserved in N % of all genes",
				"Specify N with option --consensus N; [default: 80%]",
				"Two scripts for execution in PyMol are provided:", 
				"'color_exons.py' to mark consensus exons",
				"'color_splicesites.py' to mark splice junctions of consensus exons") do |file|
				if file.start_with?("-") then 
					Helper.abort "Missing argument for --pdb <file_name>"
				end
				options[:output_format_list] << "pdb"
				vivify_hash(options, :pdb_options, :path_to_pdb, file)
				Helper.file_exist_or_die(file)
			end
			opts.on("--pdb-chain CHAIN", String,
				"Mark gene structures for chain CHAIN",
				"[Default: Use chain A]") do |opt|
				if opt.start_with?("-") then 
					Helper.abort "Missing argument for --pdb-chain <chain>"
				end
				vivify_hash(options, :pdb_options, :pdb_chain, opt)
			end
			opts.on("--pdb-ref-prot PROT", String,
				"Use protein PROT as reference for alignment with pdb sequence", 
				"[Default: First protein in alignment]") do |n|
				if n.start_with?("-") then 
					Helper.abort "Missing argument for --pdb-ref-prot <prot>"
				end
				vivify_hash(options, :pdb_options, :pdb_reference_protein, n)
			end
			opts.on("--pdb-ref-prot-struct",
				"Color only intron positions occuring in the reference protein structure.") do |opt|
				vivify_hash(options, :pdb_options, :pdb_ref_prot_struct_only, true)
			end

			opts.separator ""
			opts.on("--tree", 
				"Generate newick tree file and SVG representation") do |opt|
				options[:output_format_list] << "tree"
			end

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
			opts.on("--statistics", 
				"Output additional file with statistics about common introns",
				"To include information about taxomony, specify '--taxomony' and '--taxonomy-to-fasta' options") do 
				options[:output_format_list] << "stats"
			end

			opts.separator ""
			opts.separator "Taxonomy:"
			opts.on("--taxonomy FILE", String,
				"NCBI taxonomy database dump file FILE", 
				"OR path to extract from NCBI taxonomy:",
				"Lineage must be semicolon-separated list of taxa from root to species.") do |file|
				if file.start_with?("-") then 
					Helper.abort "Missing argument for --taxonomy <file_name>"
				end
				vivify_hash(options, :tax_options, :path_to_tax, file)
				Helper.file_exist_or_die(file)
			end
			opts.on("--taxonomy-to-fasta FILE", String,
				"Text-based file mapping gene structure file names to species names", 
				"Mandatory format:",
				"One or more genes given as semicolon-separated list and species name",
				"Delimiter between gene list and species name must be a colon",
				"The species name itself must be enclosed by double quotes like this \"SPECIES\"") do |file|
				if file.start_with?("-") then 
					Helper.abort "Missing argument for --taxonomy-to-fasta <file_name>"
				end
				vivify_hash(options, :tax_options, :path_to_tax_mapping, file)
				Helper.file_exist_or_die(file)
			end
			opts.on("--taxonomy-common-to x,y,z", Array, 
				"Mark introns common to taxa x,y,z", 
				"List must consist of at least one NCBI taxon") do |list|
				list = list.map {|str| str.capitalize}
				vivify_hash(options, :tax_options, :selected_taxa, list)
			end
			opts.on("--[no-]exclusively-in-taxa", 
				"Mark introns occuring (not) exclusively in listed taxa", 
				"Default: Not exclusively.") do |opt|
				vivify_hash(options, :tax_options, :is_exclusive, opt)
			end
			opts.on("--introns-per-taxon", 
				"Newly gained introns for every inner node in taxonomy") do |opt|
				options[:output_format_list] << "extensive_tax"
			end
			opts.on("--intron-numbers-per-taxon") do |opt|
				vivify_hash(options, :tax_options, :generate_list_intron_positios_per_taxon, true)
			end
			opts.on("--no-grep", 
				"Read the NCBI taxomony dump into RAM",
				"This will require some hundert MBs of RAM additionally",
				"Default: taxomony dump is parsed with 'grep' calls") do |opt|
				vivify_hash(options, :tax_options, :no_grep, true)
			end
			opts.on("--nice", 
				"Give grep calls to parse taxonomy dump a lower priority",
				"Please make sure to have 'nice' in your executable path when using this option") do |opt|
				vivify_hash(options, :tax_options, :nice_grep, true)
			end

			opts.separator ""
			opts.separator "Analysis and output of all or subset of data:"
			opts.on("--analyse-all-output-all", 
				"Analyse all data and provide full output [default]") do |opt|
				options[:selection][:no_selection] = true
			end
			opts.on("--analyse-all-output-selection", 
				"Analyse all data and provide text-based and graphical output for selection only", 
				"All introns are analysed, including those not present in selection") do |opt|
				options[:selection][:analyse_all_output_sel] = true
			end
			opts.on("--analyse-selection-output-selection", 
				"Analyse selected data and provide output for selection only") do |opt|
				options[:selection][:analyse_sel_output_sel] = true
			end
			opts.on("--analyse-selection-on-all-data-output-selection", 
				"Analyse intron positions of selected data in all data and provide output for selection only",
				"Introns present in selection are analysed in all data") do |opt|
				options[:selection][:analyse_sel_on_all_output_sel] = true
			end

			opts.separator("Selection criteria for data and output selection")
			opts.on("--selection-based-on-regex <\"regex\">", String, 
				"Regular expression applied on gene structure file names") do |regex|
				if regex.start_with?("-") then 
					Helper.abort "Missing argument for --selection-based-on-regex <regex>"
				end
				options[:select_by][:regex] = Regexp.new(regex)
			end
			opts.on("--selection-based-on-list x,y,z", Array, 
				"List of gene structures to be used") do |list|
				options[:select_by][:list] = list
			end
			opts.on("--selection-based-on-species <\"species\">", String, 
				"Use all gene structures associated with species", 
				"Specify also --taxonomy-to-fasta to map gene structure file names to species names") do |species|
				if species.start_with?("-") then 
					Helper.abort "Missing argument for --selection-based-on-species <species>"
				end
				options[:select_by][:species] = species
			end
			opts.on("--select-all", 
				"No selection applied (default)") do |opt|
				options[:select_by][:no_selection] = true
			end

			opts.separator ""
			opts.separator "General options:"
			opts.on("-o", "--outfile <file_name>", String, 
				"Prefix of the output file(s)", "Default: genepainter") do |file|
				if file.start_with?("-") then 
					Helper.abort "Missing argument for --outfile <file_name>"
				end
				options[:output_basename] = File.basename(file, ".*")
			end
			opts.on("--path-to-output <path>", String, 
				"Path to location for the output file(s)", 
				"Default: same location as GenePainter source files") do |path|
				if path.start_with?("-") then 
					Helper.abort "Missing argument for --path-to-output <path>"
				end
				Helper.dir_exist_or_die(path)
				Helper.dir_writable_or_die(path)
				options[:output_basepath] = path
			end

			opts.on("--range START,STOP", Array,
				"Restrict genes to range START-STOP in alignment", 
				"Might also be list if format START1,STOP1,START2,STOP2",
				"Keyword 'end' might be used to mark last position in alignment") do |range|
				# range = range.map(&:to_i)
				range = range.map do |ele| 
					if ele.downcase == "end" then 
						1.0/0.0 # Infinity value
					else
						num = ele.to_i # each non-number is mapped onto zero
						Helper.abort "Invalid argument: --range expects comma-separated list of natural numbers." if num <= 0
						Helper.human2ruby_counting( ele.to_i ) 
					end
				end

				# sanity check
				# range should be 2 (or multiple of two) natural numbers
				if range.size % 2 != 0 then
					Helper.abort "Invalid argument: --range expects two (or multiple of two) comma-separated natural numbers without spaces"
				end
				# range_start should be smaller than r_stop
				range.each_slice(2) do |r_start, r_stop|
					if r_stop != "end" && r_start >= r_stop then 
						Helper.abort "Invalid range definition in #{range}: start #{r_start} must be less than stop #{r_stop}"
					end
				end
				# _end sanity check
			options[:range][:reverse_position] = range.sort.reverse
				# options[:range][:position] = {
					# position: range, is_delete_range: false
					# }			
			end

			opts.on("--[no-]delete-range", 
				"(Not) Delete specified range") do |opt|
				options[:range][:is_delete_range] = opt
				# vivify_hash(options, :range, :is_delete_range, opt)
			end

			opts.on("--keep-common-gaps", 
				"Keep common gaps in alignment") do
				options[:ignore_common_gaps] = false
			end

			opts.on("--no-best-position-introns", 
				"Plot introns always onto beginning of a gap",
				"Default: Align introns if their position differs by alignment gaps only") do 
				options[:best_intron_pos] = false
			end

			opts.on("--[no-]separate-introns-in-textbased-output", 
				"(Not) Separate each consecutive pair of introns by an exon placeholder in text-based output formats.",
				"Default: Separate introns unless the output lines get too long.") do |opt|
				options[:is_long_text_based_output] = opt
			end
	
			opts.separator ""
			opts.on_tail("-h", "--help", "Show this message") do 
				# delete all hidden params from help
				show_help(opt_parser.help)
				exit
			end

		end # optionparser

		## main part of function ##
		if args.empty? then
			# display help and exit if program is called without any argument
			show_help(opt_parser.help) 
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

		# genereate the path to output from prefix name and path
		if options[:output_basepath] != "" then 
			options[:path_to_output] = File.join(options[:output_basepath], options[:output_basename])
		else
			options[:path_to_output] = options[:output_basename]
		end
		options[:path_to_log] = options[:path_to_output] + ".log"
		options.delete(:output_basename)
		options.delete(:output_basepath)

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
		# if taxonomy should be considered, need taxonomy dump file, mapping fasta to tax file and list of taxa
		if options[:tax_options].any? ||
			options[:output_format_list].include?("extensive_tax") || 
			options[:output_format_list].include?("tree") then 
			if ! (
				options[:tax_options][:path_to_tax] &&
				options[:tax_options][:path_to_tax_mapping]	) then 

				Helper.abort "Mandatory argument for taxonomy is missing: Specify --taxonomy FILE and --taxonomy-to-fasta FILE"
			end
		end

		# dont make taxonomy mandatory for statistiscs option!
		if options[:output_format_list].include?("stats") && 
			(options[:tax_options][:path_to_tax] || options[:tax_options][:path_to_tax_mapping]) then 
			if ! ( options[:tax_options][:path_to_tax] && options[:tax_options][:path_to_tax_mapping] ) then
				Helper.abort "Mandatory argument for --statistics in combination with taxonomy is missing: Specify --taxonomy FILE and --taxonomy-to-fasta FILE"
			end
		end

		# only one selection type can be selected
		# need to be combined with one selection criterium
		if options[:selection][:no_selection] && 
			(options[:select_by][:regex] || options[:select_by][:list] || options[:select_by][:species]) then 
				Helper.abort "Invalid usage of --analyse-all-output-all. Cannot be combined with selection criteria."
		end
		if options[:select_by][:no_selection] && 
			(options[:selection][:analyse_sel_output_sel] || 
				options[:selection][:analyse_sel_on_all_output_sel] || 
				options[:selection][:analyse_all_output_sel]) then 
			Helper.abort "Invalid usage of --select-all. Cannot be combined with analysis or output of subset of data."
		end
		n_specified_selections = 0
		options[:selection].each do |opt, val|
			n_specified_selections += 1 if val
			Helper.abort "Invalid usage of analysis or output selection options. Specify not more than one criterium." if n_specified_selections > 2
		end
		n_specified_selections -= 1 if options[:selection].delete(:no_selection)
		n_specified_criteria = 0
		options[:select_by].each do |opt, val|
			n_specified_criteria += 1 if val
			Helper.abort "Invalid usage of selection criteria. Specify not more than one criterium." if n_specified_criteria > 2
		end
		n_specified_criteria -= 1 if options[:select_by][:no_selection]
		if n_specified_selections == 1 && n_specified_criteria != 1 then 
			Helper.abort "Invalid usage of analysis or output selection options. Must be combined with selection criterium."
		end
		if n_specified_criteria == 1 && n_specified_selections != 1 then 
			Helper.abort "Invalid usage of selection criterium. Must be combined with analysis or output selection options."
		end
		if options[:select_by][:species] && ! options[:tax_options][:path_to_tax_mapping] then 
			Helper.abort "Invalid usage of --selection-based-on-species: Must be combined with --taxonomy-to-fasta"
		end

		return options

		# use the own format of fatal error messages!				
		rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption, OptionParser::AmbiguousOption => exc
			Helper.abort exc
	end # parse()

	def self.vivify_hash(hash, outer_key, inner_key, inner_val)
		if ! hash[outer_key] then 
			hash[outer_key] = {}
		end
		hash[outer_key][inner_key] = inner_val
	end

	def self.show_help(opts)
		opts_array = opts.to_s.split("\n")
		@hidden_params.each do |hidden_param|
			opts_array.delete_if{ |line| line =~ /#{hidden_param}/ }
		end
		puts opts_array.join("\n")
	end	

end
