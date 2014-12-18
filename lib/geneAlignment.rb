# aligns gene objecs

class GeneAlignment
	attr_reader :exon_placeholder, :intron_placeholder

	## output parameter definition 
	@consensus_val = "1.0"
	def self.consensus_val=(val)
		@consensus_val = val.to_s
	end
	def self.consensus_val
		@consensus_val
	end
	def self.max_length_gene_name
		return 20
	end
	def self.max_introns_for_long_reduced_output
		return 50
	end
	def self.suffix_structure_in_alignment
		return "_structure"
	end
	def self.merged_structure_name
		return "Merged"
	end
	def self.consensus_structure_name
		return "Consensus_#{consensus_val}"
	end
	def self.taxonomy_structure_name
		return "Taxonomy"
	end
	## end output parameter definition

	# input
	# Array genes: gene objects, must have an aligned sequence
	# Float or false: calculate a consensus pattern (conserved in % val sequences)
	# Boolean: calculate a merged pattern
	# Boolean: calculate statistics
	# Hash taxonomy_options: genes_within_taxa: Array [subset of gene.names (belong to genes objects)], is_exclusive: Boolean [introns exclusive for selected taxa]
	# Boolean: reduced exon-intron pattern contains no common gaps at all or one of each series
	# Boolean: output full instead of reduced exon-intron pattern
	def initialize(genes, consensus_val, is_merged_pattern, taxonomy_options, sep_introns_in_plaintext_output, use_full_patterns_in_plaintext_output)

		@genes = genes # containing all genes, independent of genestructures

		@exon_placeholder = "-"
		@intron_placeholder = nil # defaults to intron phase

		@is_separate_introns_in_textbased_output = set_separate_introns_value(sep_introns_in_plaintext_output)
		@is_use_full_pattern_in_plaintext_output = use_full_patterns_in_plaintext_output

		# overwritten by method convert_to_exon_intron_pattern
		@ind_consensus_pattern = nil
		@ind_merged_pattern = nil
		@ind_tax_pattern = nil

		@aligned_genestructures, @names_aligned_genestructures, @stats_per_intron_pos = convert_to_exon_intron_pattern 
		@additional_structures = convert_to_special_exon_intron_pattern(consensus_val, is_merged_pattern, taxonomy_options)
		@reduced_aligned_genestructures, @reduced_additional_structures = reduce_exon_intron_pattern() # reduce gene structures to "needed" parts 

		@n_structures = calc_n_structures
	rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "initialize gene alignment", exp
    	Helper.abort "Cannot create gene alignment."
	end

	def calc_n_structures
		@aligned_genestructures.size + @additional_structures.size
	end

	def set_separate_introns_value(sep_introns_in_plaintext_output)
		boolean = sep_introns_in_plaintext_output 
		if boolean.nil? then 
			# if not explicitly set, set value in relation to number of intron positions (to keep output small)
			if @stats_per_intron_pos.keys.size > self.class.max_introns_for_long_reduced_output then 
				boolean = false
			end	
		end
		return boolean
	end

	def get_gene_obj_by_name(searched_name)
		return @genes.find{ |g| g.name == searched_name }
	end

	# align genestructures by plotting them onto the multiple sequence alignment
	# exon representation: "-", intron representation: phase
	# output
	# Array of strings: exon intron pattern plotted onto the aligned sequence
	# Array of strings: gene names in same order as exon intron patterns
	# Hash: information about each intron position found
	def convert_to_exon_intron_pattern
		# collect the exon_intron_patterns of each gene, also for merged/consensus

		patterns = Array.new(@genes.size)
		names = Array.new(@genes.size)

		# number of introns, phase, and genes it occurs in for each intron position
		# taxon_first_found: the last common ancestor of all species, genes are associated with
		intron_pos_with_annotation = Hash.new { |h,k| h[k] = { n_introns: 0, phase: "?", genes: [], taxon_first_found: "" } }

		@genes.each_with_index do |gene, ind|
			patterns[ind] = gene.plot_intron_phases_onto_aligned_seq(@exon_placeholder, @intron_placeholder)
			names[ind] = gene.name

			update_annotation_per_intron_position(gene, intron_pos_with_annotation)
		end

		# attach the first taxon the intron is found in to each intron position
		add_taxon_first_found_to_annotation_per_intron_position(intron_pos_with_annotation)

		return patterns, names, intron_pos_with_annotation
	end
	def convert_to_special_exon_intron_pattern(consensus_val, is_merged_pattern, taxonomy_options)
		# generate three add. patterns, delete unused ones
		patterns = []

		pattern_length = @aligned_genestructures[0].size

		# parse taxonomy input
		genes_belonging_to_selected_taxa = taxonomy_options[:genes_belonging_to_selected_taxa]
		is_intron_exclusive_for_taxa = taxonomy_options[:is_intron_exclusive_for_selected_taxa]

		if consensus_val then 
			min_n_introns_per_pos = consensus_val * @genes.size
			patterns << generate_consensus_profile( @stats_per_intron_pos, min_n_introns_per_pos , pattern_length )
			@ind_consensus_pattern = patterns.size - 1 # indices start with 0
			self.class.consensus_val = consensus_val # set class variable needed for export functions
		end
		if is_merged_pattern then
			patterns << generate_merged_profile( @stats_per_intron_pos, pattern_length )
			@ind_merged_pattern = patterns.size - 1 # indices start with 0
		end
		if genes_belonging_to_selected_taxa.any? then
			patterns << generate_taxonomy_profile( @stats_per_intron_pos, 
				genes_belonging_to_selected_taxa, is_intron_exclusive_for_taxa, 
				pattern_length )
			@ind_tax_pattern = patterns.size - 1 # indices start with 0
		end
		return patterns
	end

	# update counts, phase and occurence-list per intron position
	def update_annotation_per_intron_position(gene, intron_pos_with_annotation)

		gene.get_all_introns_with_phase.each do |pos, phase|
			pos_phase = Intron.merge_position_and_phase(pos, phase)

			intron_pos_with_annotation[pos_phase][:phase] = phase # yes, this is redundant, but also convenient :-)
			intron_pos_with_annotation[pos_phase][:n_introns] += 1
			intron_pos_with_annotation[pos_phase][:genes].push gene.name

		end
	end

	def add_taxon_first_found_to_annotation_per_intron_position(intron_pos_with_annotation)
		intron_pos_with_annotation.each do |pos, data|
			genes_having_this_intron = data[:genes]
			lineage_common_to_all = []

			# collect the lineage that is common to all
			genes_having_this_intron.each do |gene_name|
				gene = get_gene_obj_by_name(gene_name) # independent from output-selection of genes

				if gene.has_taxonomic_information then 
					# taxonomic lineage of gene is known
					if lineage_common_to_all.empty? then 
						# no lineage_common_to_all found so far, set it to lineage of this gene
						lineage_common_to_all = gene.get_lineage_root_to_first_uniq_ancestor_of_species
					else
						# find parts of lineage_common_to_all that also occur in this gene
						lineage_common_to_all = gene.common_lineage_with_other_gene(lineage_common_to_all)
					end
				else
					# taxonomic lineage is unknown
					# nothing to do.
				end
			end
			# taxon_first_found is nil if no lineage is known
			intron_pos_with_annotation[pos][:taxon_first_found] = lineage_common_to_all.last # lineage starts with root
		end
	end


	# methods for "statistics per intron position and statistics: merged/consensus/taxonomy profile"

	def generate_consensus_profile(intron_pos_with_annotation, min_n_introns, pattern_length)
		consensus_pattern = get_empty_pattern(pattern_length)

		intron_pos_with_annotation.each do |intronpos, info|
			# check if the intron occurs often enough
			if info[:n_introns] >= min_n_introns then 
				consensus_pattern[intronpos] = info[:phase]
			end
		end
		return consensus_pattern
	end
	def generate_merged_profile(intron_pos_with_annotation, pattern_length)
		merged_pattern = get_empty_pattern(pattern_length)

		intron_pos_with_annotation.each do |intronpos, info|
			merged_pattern[intronpos] = info[:phase]
		end
		return merged_pattern
	end
	def generate_taxonomy_profile( intron_pos_with_annotation, genes_in_selected_taxa, is_intron_exclusive_for_taxa, pattern_length )
		tax_pattern = get_empty_pattern(pattern_length)

		intron_pos_with_annotation.each do |intronpos, info|
			if ( is_intron_exclusive_for_taxa && info[:genes].is_subset?(genes_in_selected_taxa) ) ||
				( ! is_intron_exclusive_for_taxa && info[:genes].is_overlapping_set?(genes_in_selected_taxa) ) then 
				tax_pattern[intronpos] = info[:phase]
			end
		end
		return tax_pattern
	end
	def get_empty_pattern(len)
		return @exon_placeholder * len
	end
	# end methods for statistics

	# delete genes and gene structures not selected for output
	# optional: delete also all introns occuring in genes not selected for output
	# generates reduced_aligned_genestructres new to really remove all uncessary exon-placeholders
	def reduce_gene_set_for_output(selected_for_output, is_delete_introns_not_occuring_in_selection, consensus_val, is_merged_pattern, taxonomy_options)
		n_selected_genes = selected_for_output.size

		new_names = Array.new(n_selected_genes)
		new_aligned_genestructures = Array.new(n_selected_genes)

		# generate reduced set of gene names, gene structures
		selected_for_output.each_with_index do |selected_gene_name, ind|
			# find index of this gene in old data set
			ind_old = @names_aligned_genestructures.index(selected_gene_name)

			new_names[ind] = selected_gene_name
			new_aligned_genestructures[ind] = @aligned_genestructures[ind_old]
		end

		# generate reduced set of intron positions (if neccessary)
		if is_delete_introns_not_occuring_in_selection then 
			@stats_per_intron_pos.delete_if do |intron_pos, intron_info|
				intron_info[:genes].is_disjoint_set?(selected_for_output)
			end

			@is_reduced_intronpositions = true
		end

		# overwrite old class variables
		@names_aligned_genestructures = new_names
		@aligned_genestructures = new_aligned_genestructures
		# generate special patterns, automatically overwrites @ind_ -variables
		@additional_structures = convert_to_special_exon_intron_pattern(consensus_val, is_merged_pattern, taxonomy_options)
		@reduced_aligned_genestructures, @reduced_additional_structures = reduce_exon_intron_pattern()

		@n_structures = calc_n_structures
    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "reduce_gene_set_for_output", exp
    	Helper.abort "Cannot reduce gene set for output."
	end	

	def reduce_exon_intron_pattern(opts={})
		is_del_all_common_gaps = opts[:del_all_common_gaps] || false
		is_sep_introns = opts[:sep_introns_in_plaintext_output] || @is_separate_introns_in_textbased_output
		input_array = opts[:input] || [ @aligned_genestructures,@additional_structures ]

		all_patterns = Sequence.remove_common_gaps(
			input_array.flatten,
			{delete_all_common_gaps: is_del_all_common_gaps,
			ensure_common_gap_between_consecutive_non_gaps: is_sep_introns,
			gap_symbol: @exon_placeholder} # use placeholder as gap-symbol
			)

		last_ind_gene_pattern = @aligned_genestructures.size - 1
		gene_patterns = all_patterns[0..last_ind_gene_pattern]
		add_patterns = all_patterns[(last_ind_gene_pattern+1)..-1]
		return gene_patterns, add_patterns
	end

	def convert_string_to_chopped_fasta_header(str)
		# add ">" at beginning of string
		# left justify string to certain length
		(">" << str)[0...self.class.max_length_gene_name].ljust(self.class.max_length_gene_name)
	end

	def replace_exon_intron_placeholders_in_structure(struct, exon_placeholder_final, intron_placeholder_final)

		# replace default placeholders for exons and introns by the requested ones
		# order does matter: start with introns, in case of binary output
		# (otherwise: all exons would 0, afterwards all 0,1,2 become 1)
		if intron_placeholder_final != @intron_placeholder then
			struct = struct.gsub(intron_placeholder_regexp, intron_placeholder_final)
		end

		if exon_placeholder_final != @exon_placeholder then
			struct = struct.gsub(@exon_placeholder, exon_placeholder_final)
		end

		return struct
	end
	def intron_placeholder_regexp
		# this is neccessary, as the default placeholder is "nil"
		if @intron_placeholder then
			return Regexp.new( Regexp.escape( @intron_placeholder ) )
		else
			return Regexp.new( "[0|1|2|?]" )
		end
	end

	def export_as_alignment_with_introns

		# output contains sequence and structure for each gene, but no sequence for merged/conserved structure
		output = Array.new(@genes.size + @n_structures)

		n_additional_structs = 0
		@additional_structures.each_with_index do |struct, ind|
			if ind == @ind_merged_pattern then 
				name_struct = self.class.merged_structure_name
			elsif ind == @ind_consensus_pattern
				name_struct = self.class.consensus_structure_name
			elsif ind == @ind_tax_pattern
				name_struct = self.class.taxonomy_structure_name
			end
			output[ind] = Sequence.convert_strings_to_fasta(name_struct, struct)
			n_additional_structs += 1
		end

		@aligned_genestructures.each_with_index do |struct, ind|
			# its the structure of a @aligned_gene

			index_output_struct = ind * 2 + 1 + n_additional_structs # add offset for additional structures
			index_output_seq = ind * 2 + n_additional_structs
			gene_name = @names_aligned_genestructures[ind]
			gene = get_gene_obj_by_name(gene_name)
			output[index_output_seq] = Sequence.convert_strings_to_fasta(gene_name, gene.aligned_seq)
			name_struct = gene_name + self.class.suffix_structure_in_alignment

			output[index_output_struct] = Sequence.convert_strings_to_fasta(name_struct, struct)
		end

		return output.join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_alignment_with_introns", exp
    	throw(:error)
	end

	def export_as_binary_alignment

		# output contains binary gene structure, in fasta format
		output = Array.new(@n_structures)

		exon_placeholder_output = "0"
		intron_placeholder_output = "1"

		# output should not contain any columns containing zeros only
		red_structs, red_add_structs = reduce_exon_intron_pattern({del_all_common_gaps: true}) # true: is delete all common gaps

		n_additional_structs = 0
		red_add_structs.each_with_index do |struct, ind|
			if ind == @ind_merged_pattern then 
				name = self.class.merged_structure_name
			elsif ind == @ind_consensus_pattern
				name = self.class.consensus_structure_name
			elsif ind == @ind_tax_pattern
				name = self.class.taxonomy_structure_name
			end
			name = convert_string_to_chopped_fasta_header ( name )
			struct = replace_exon_intron_placeholders_in_structure(struct, exon_placeholder_output, intron_placeholder_output)

			output[ind] = Sequence.convert_strings_to_fasta( name, struct )	
			n_additional_structs += 1
		end

		red_structs.each_with_index do |struct, ind|
			# its a gene
			name = @names_aligned_genestructures[ind]
			struct = replace_exon_intron_placeholders_in_structure(struct, exon_placeholder_output, intron_placeholder_output)

			output[ind + n_additional_structs] = Sequence.convert_strings_to_fasta( name, struct )		
		end
		return output.join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_binary_alignment", exp
    	throw(:error)
	end

	def export_as_plain_txt(exon_placeholder_output, intron_placeholder_output )
	
		output = Array.new(@n_structures)
		n_additional_structs = 0
		
		# containing either reduced or full structures
		structs_for_output = [] 
		additional_structs_for_output = []

		if @is_use_full_pattern_in_plaintext_output then 
			structs_for_output = @aligned_genestructures
			additional_structs_for_output = @additional_structures
		else
			structs_for_output = @reduced_aligned_genestructures
			additional_structs_for_output = @reduced_additional_structures
		end

		additional_structs_for_output.each_with_index do |struct, ind|
			if ind == @ind_merged_pattern then 
				name = self.class.merged_structure_name
			elsif ind == @ind_consensus_pattern
				name = self.class.consensus_structure_name
			elsif ind == @ind_tax_pattern
				name = self.class.taxonomy_structure_name
			end

			name = convert_string_to_chopped_fasta_header ( name )
			struct = replace_exon_intron_placeholders_in_structure(struct, exon_placeholder_output, intron_placeholder_output)
			
			output[ind] = [name, struct].join("")	
			n_additional_structs += 1
		end

		structs_for_output.each_with_index do |struct, ind|
			# its a gene
			name = @names_aligned_genestructures[ind]

			name = convert_string_to_chopped_fasta_header ( name )
			struct = replace_exon_intron_placeholders_in_structure( struct, exon_placeholder_output, intron_placeholder_output )

			output[ind+n_additional_structs] = [name, struct].join("")
					
		end
		return output.join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_plain_txt", exp
    	throw(:error)
	end

	# generates three-fold output
	# 1) original pattern with all fuzzy introns (within window_size) marked
	# 2) list of all positions that will be merged
	# 3) fuzzy pattern with all fuzzy introns merged to one position
	def export_as_plain_text_with_fuzzy_intron_pos(exon_placeholder_output, intron_placeholder_output, window_size)

		sorted_intron_positions = @stats_per_intron_pos.keys.sort
		matched_fuzzy_positions = {} # key: original position, value: mapped position
		fuzzy_intron_placeholder_output = "*"

		# init 'results' variables
		output = [
			["Original exon-intron pattern, fuzzy introns (window size #{window_size}) marked by #{fuzzy_intron_placeholder_output}"],
			["\n", "List of fuzzy introns"],
			["\n", "Exon-intron pattern with fuzzy introns written at same place"]
		] # contains unfuzzy pattern with fuzzy introns marked, fuzzy_pattern and info
		unfuzzy_pattern_with_fuzzy_marked = @aligned_genestructures.map do |ele| 
			replace_exon_intron_placeholders_in_structure(
				ele.dup,
				exon_placeholder_output,
				intron_placeholder_output
			)
		end
		fuzzy_pattern = unfuzzy_pattern_with_fuzzy_marked.map {|ele| ele.dup}
		info_per_intronpos = [ "Intron position\tMapped intron positions" ]

		sorted_intron_positions.each_cons(2) do |pos1, pos2|

			if Intron.are_introns_within_range(pos1, pos2, window_size) then
				# introns are fuzzy: their position is within window_size
				# map them to first 'real' position

				# intronpositions that occur in the same gene are not fuzzy!
				genes_common_to_pos1_and_pos2 = @stats_per_intron_pos[pos1][:genes] & @stats_per_intron_pos[pos2][:genes]

				if map_intron_onto_pos = matched_fuzzy_positions[pos1] then 
					matched_fuzzy_positions[pos1] = map_intron_onto_pos
					matched_fuzzy_positions[pos2] = map_intron_onto_pos
				else
					map_intron_onto_pos = pos1
					matched_fuzzy_positions[pos2] = pos1
				end

				(@stats_per_intron_pos[pos1][:genes] - genes_common_to_pos1_and_pos2).each do |gene_name|
					ind = @names_aligned_genestructures.index(gene_name)
					next if ! ind # gene does not belong to output

					unfuzzy_pattern_with_fuzzy_marked[ind][pos1] = fuzzy_intron_placeholder_output
					if map_intron_onto_pos != pos1 then 
						fuzzy_pattern[ind][pos1] = exon_placeholder_output
						fuzzy_pattern[ind][map_intron_onto_pos] = intron_placeholder_output
					end
				end

				(@stats_per_intron_pos[pos2][:genes] - genes_common_to_pos1_and_pos2).each do |gene_name|
					ind = @names_aligned_genestructures.index(gene_name)
					next if ! ind # gene does not belong to output

					unfuzzy_pattern_with_fuzzy_marked[ind][pos2] = fuzzy_intron_placeholder_output
					fuzzy_pattern[ind][pos2] = exon_placeholder_output
					fuzzy_pattern[ind][map_intron_onto_pos] = intron_placeholder_output
				end
			else
				# no fuzzy intron
			end

		end

		# reduce generated patterns
		reduced_unfuzzy_pattern_with_legend = exon_intron_pattern_with_merged_pattern_spaces_and_intron_numbers(
			unfuzzy_pattern_with_fuzzy_marked, exon_placeholder_output, intron_placeholder_output)
		reduced_fuzzy_pattern = reduce_exon_intron_pattern( {input: fuzzy_pattern} ).first # first output is pattern itself

		# add gene names to fuzzy pattern
		reduced_fuzzy_pattern.each_with_index do |struct, ind|
			name = @names_aligned_genestructures[ind]

			name = convert_string_to_chopped_fasta_header ( name )
			reduced_fuzzy_pattern[ind] = [name, struct].join("")
		end

		# generate info which introns are fuzzy
		matched_fuzzy_positions.values.sort.uniq.each do |mapped_pos|
			original_positions = matched_fuzzy_positions.select{|k,v| v==mapped_pos}.keys

			ind_mapped_intron_pos = Helper.ruby2human_counting( sorted_intron_positions.index(mapped_pos) ) # start count with 1 !
			ind_original_positions = original_positions.collect {|pos| Helper.ruby2human_counting( sorted_intron_positions.index(pos) ) }
			info_per_intronpos.push( "#{ind_mapped_intron_pos}\t#{ind_original_positions.join(", ")}" )
		end

		# add legend with intron numbers to unfuzzy pattern and this pattern to output
		output[0] << reduced_unfuzzy_pattern_with_legend.join("\n")
		output[1] << info_per_intronpos.join("\n")
		output[2] << reduced_fuzzy_pattern.join("\n")

		return output.join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_plain_text_with_fuzzy_intron_pos", exp
    	throw(:error)
	end

	def export_as_svg(options)
		
		gene_objs_for_output = @genes.select do |gene|
			@names_aligned_genestructures.include?(gene.name)
		end

		# image can be in format "normal" or "reduced"
		# also, both formats can be requested at the same time
		output_normal_format = nil
		output_reduced_format = nil
		is_nested_svg_elements = options[:generate_nested_svg] # creates nested svg elements for easier positioning

		if options[:reduced] || options[:both] then 
			is_default_output = false

			# prepare data
			genealignment2svg_obj = GeneAlignment2svg.new(gene_objs_for_output, is_default_output, is_nested_svg_elements)

			# draw genes
			output_reduced_format = genealignment2svg_obj.create_svg
		end
		if ! options[:reduced] || options[:both] then 
			is_default_output = true

			# prepare data
			genealignment2svg_obj = GeneAlignment2svg.new(gene_objs_for_output, is_default_output, is_nested_svg_elements)

			# draw genes
			output_normal_format = genealignment2svg_obj.create_svg
		end

		return output_normal_format, output_reduced_format

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_svg", exp
    	throw(:error)
	end

	def export_as_svg_only_merged_pattern(options)

		gene_objs_for_output = @genes.select do |gene|
			@names_aligned_genestructures.include?(gene.name)
		end
		
		# image can be in format "normal" or "reduced"
		# also, both formats can be requested at the same time
		output_normal_format = nil
		output_reduced_format = nil
		is_nested_svg_elements = options[:generate_nested_svg] # creates nested svg elements for easier positioning

		if options[:reduced] || options[:both] then 
			is_default_output = false

			# prepare data
			genealignment2svg_obj = GeneAlignment2svg.new(gene_objs_for_output, is_default_output, is_nested_svg_elements)

			# draw merged gene
			output_reduced_format = genealignment2svg_obj.create_svg_merged_genestructure

		end
		if ! options[:reduced] || options[:both] then 
			is_default_output = true

			# prepare data
			genealignment2svg_obj = GeneAlignment2svg.new(gene_objs_for_output, is_default_output, is_nested_svg_elements)

			# draw merged gene
			output_normal_format = genealignment2svg_obj.create_svg_merged_genestructure
		end
		return output_normal_format, output_reduced_format

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_svg", exp
    	throw(:error)
	end

	def export_as_pdb(options)
		output = []
		Helper.log "\nMap gene structures onto PDB #{options[:path_to_pdb]}"
		GeneAlignment2pdb.log_howtouse_pythonscripts

		# prepare data
		ref_seq = "" # the aligned reference sequence
		ref_struct = "" # the genestructure plotted onto the ref_seq
		ref_name = "" # name of selected reference gene

		# default: use first seq in alignment as reference
		# if someother gene is specified via command line, overwrite sequence index
		ind_ref_name = 0
		if options[:pdb_reference_protein] then 
			# use specified seq
			name = options[:pdb_reference_protein]
			if name.start_with?(">") then
				name = name[1..-1]
			end
			ind_ref_name = @names_aligned_genestructures.index(name)
			if ! ind_ref_name then
				Helper.log "Cannot match gene #{options[:pdb_reference_protein]}."
				Helper.log "Use the first gene in alignment instead."		
				ind_ref_name = 0
			end
		end

		ref_name = @names_aligned_genestructures[ind_ref_name]
		ref_seq = get_gene_obj_by_name(ref_name).aligned_seq

		# get gene structure
		ind_ref_struct = nil
		if options[:pdb_ref_prot_struct_only]
			# use structure of reference sequnce
			ref_struct = @aligned_genestructures[ind_ref_name]
		else
			# use one of the additional patterns?
			if @ind_merged_pattern && ind_ref_struct.nil? then 
				# use the merged pattern
				ind_ref_struct = @ind_merged_pattern
				Helper.log "Selected merged pattern for mapping onto PDB"
			end
			if @ind_consensus_pattern && ind_ref_struct.nil? then 
				ind_ref_struct = @ind_consensus_pattern
				Helper.log "Selected consensus pattern for mapping onto PDB"
			end
			if @ind_tax_pattern && ind_ref_struct.nil? then 
				ind_ref_struct = @ind_tax_pattern
				Helper.log "Selected taxonomic pattern for mapping onto PDB"
			end
			ref_struct = @additional_structures[ind_ref_struct]
		end
		if ref_struct.empty? then 
			# default: use first structure 
			ref_struct = @aligned_genestructures[0]
		end

		genealignment2pdb_obj = GeneAlignment2pdb.new( ref_seq, ref_struct, options )

		# plot gene structures onto pdb
		output1, output2 = genealignment2pdb_obj.map_genestructure_onto_pdb

		genealignment2pdb_obj.log_alignment(ref_name, options[:path_to_pdb])

		return output1, output2

    # rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    # 	Helper.log_error "export_as_pdb", exp
    # 	throw(:error)		
	end

	def export_as_tree

		# init vars for annotating the tree
		taxa_with_alternative_names = {} # final annotation
		all_taxa_with_statistics = {} # statistics about # species, # genes, # intronpositions (for each leaf); literally for all taxa
		taxa_with_gain_loss = {} # statistics about # gained introns, # lost introns; only for those taxa where gain/loss happend

		all_lineages = []
		common_ancestors = nil

		# collect all lineages found and statistics like the number of genes and intronpositions per taxon
		@genes.each do |gene|
			if gene.has_taxonomic_information then 
				
				all_lineages |= [gene.taxonomic_lineage]
				species = gene.taxonomic_lineage.last

				if common_ancestors then 
					common_ancestors = gene.taxonomic_lineage & common_ancestors
				else
					common_ancestors = gene.taxonomic_lineage
				end

				# collect species, gene and intron position counts for each leaf
				gene.get_lineage_root_to_first_uniq_ancestor_of_species.each do |taxon|
					update_statistics(all_taxa_with_statistics, taxon, gene)
				end
			end
		end

		taxa_intron_first_found = []
		# collec all taxa, at which an intron is first found and add update gain-counts
		@stats_per_intron_pos.each do |intronpos, introninfo|
			if taxon = introninfo[:taxon_first_found] then 
				taxa_intron_first_found |= [taxon]
			end
		end

		# stop execution here if there is no taxonomic information !!!
		throw :no_taxonomy if all_lineages.empty? || taxa_intron_first_found.empty? || common_ancestors.nil? 

		last_common_ancestor = find_last_common_ancestor(all_lineages, common_ancestors, taxa_intron_first_found)

		tree_obj = Tree.new(all_lineages, last_common_ancestor, taxa_intron_first_found)

		# map genes onto leaves
		genes_with_corresponding_leaves = map_genes_onto_leaves(tree_obj.get_leaves)

		# collect gain/loss counts
		@stats_per_intron_pos.each do |intronpos, introninfo|
			if taxon = introninfo[:taxon_first_found] then 
				# genes correspond to which leaves???
				leaves_with_intron = introninfo[:genes].collect{ |gene| genes_with_corresponding_leaves[gene] }.compact.uniq

				nodes_without_intron = tree_obj.find_nodes_without_intronposition(taxon, leaves_with_intron)

				update_gain(taxa_with_gain_loss, taxon)
				update_losses(taxa_with_gain_loss, nodes_without_intron)
			end
		end

		# use collected data as alternative names
		# format:
		# name[_gain/loss]?[_stats]?
		# special features:
		# - nodes without any gain or loss get 0/0 counts
		# - nodes with 0 gains due to the removal of intron positions: 0 gain is replaced by "-"
		tree_obj.get_all_nodes.each do |node|
			if tree_obj.is_leaf?(node) && all_taxa_with_statistics[node] then 
				taxa_with_alternative_names[node] = prettify_statisics( all_taxa_with_statistics, node)
			end

			if taxa_with_gain_loss[node] then 
				# introns gained/lost at this node
				data = taxa_with_gain_loss[node]
				if @is_reduced_intronpositions then 
					# intron positions were removed, so 0 gains is incorrect: they were just not determinded. do not display them
					is_display_zero_gain = false 
				else
					is_display_zero_gain = true # display 0 gains as this is the correct information here
				end
			else
				# nothing happend, name it like this!
				data = init_gain_loss
				is_display_zero_gain = true # display 0 gains as this is the correct information here
			end

			# update alternative-names list
			if taxa_with_alternative_names[node] then 
				taxa_with_alternative_names[node] += "_" + prettify_gain_loss_without_taxonname( data, is_display_zero_gain )
			else
				taxa_with_alternative_names[node] = prettify_gain_loss( data, node, is_display_zero_gain )
			end
		end
		output = tree_obj.export_tree(taxa_with_alternative_names)
		return output

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_tree", exp
    	throw(:error)	
	end
	# --- begin helper methods for export_as_tree
	def find_last_common_ancestor(all_lineages, common_ancestors, taxa_intron_first_found)
		# # 1) simple case: use the last of all common ancestors
		# # 2) more difficult case: all needed taxa belong to one lineage and there is only this one lineage -> use the first needed taxon as root

		if all_lineages.size == 1 && taxa_intron_first_found.is_subset?(all_lineages[0]) then 
			# case 2)
			return all_lineages[0].intersection(taxa_intron_first_found).first # the first (nearest to root) taxon really needed
		else
			# case 1)
			return common_ancestors.last # the last (farest to root) taxon common to all
		end		
	end

	def map_genes_onto_leaves(all_leaves)
		genes_with_leaves = {}
		@genes.each do |gene|
			if gene.has_taxonomic_information then 
				gene.taxonomic_lineage.reverse.each do |taxon|
					if all_leaves.include?(taxon) then 
						# found the corresponding leaf to this gene
						genes_with_leaves[gene.name] = taxon
						break # continue with next gene
					end
				end
			end # gene.has_tax_info
		end # @genes.each
		return genes_with_leaves
	end
	def update_statistics(statistics, key, gene)
		if !statistics[key] then 
			# init value 
			statistics[key] = {species: [], genes: [], intronpos: []}
		end
		statistics[key][:species] |= [ gene.taxonomic_lineage.last ]
		statistics[key][:genes] |= [ gene.name] 
		statistics[key][:intronpos] |= gene.get_all_intronpositions_merged_with_phase
	end
	def prettify_statisics(statistics, key)
		data = statistics[key]
		Helper.sanitize_taxon_name(key) + ".#" + data[:species].size.to_s + ".#" + data[:genes].size.to_s + ".#" + data[:intronpos].size.to_s
	end
	def init_gain_loss
		{gain: 0, loss: 0}
	end
	def update_gain(statistics, key)
		if ! statistics[key] then 
			statistics[key] = init_gain_loss
		end
		statistics[key][:gain] += 1
	end
	def update_losses(statistics, keys)
		keys.each do |key|
			if ! statistics[key] then 
				statistics[key] = init_gain_loss
			end
			statistics[key][:loss] += 1
		end
	end
	def prettify_gain_loss(data, key, is_display_zero_gain)
		Helper.sanitize_taxon_name(key) + "_" + prettify_gain_loss_without_taxonname(data, is_display_zero_gain)
	end
	def prettify_gain_loss_without_taxonname(data, is_display_zero_gain)
		if ( ! is_display_zero_gain)  && data[:gain] == 0 then 
			gain = "-"
		else
			gain = data[:gain].to_s
		end
		gain + "green_" + data[:loss].to_s + "red"
	end
	# --- end helper methods for export_as_tree


	# creates imaginative patterns for each taxon at which an intron is found first
	# this imaginative patterns contain all introns that can be assigned to that taxon
	def export_as_taxonomy
		taxon_first_found_with_pattern = {}
		pattern_length = @aligned_genestructures.first.size

		@stats_per_intron_pos.each do |intronpos, introninfo|
			if taxon = introninfo[:taxon_first_found] then 
				if ! taxon_first_found_with_pattern[taxon] then 
					taxon_first_found_with_pattern[taxon] = get_empty_pattern(pattern_length)
				end

				# add this intronposition to pattern
				taxon_first_found_with_pattern[taxon][intronpos] = introninfo[:phase]
			end
		end

		# stop execution here if there is no taxonomic information!!!
		throw :no_taxonomy if taxon_first_found_with_pattern.empty?

		output = Array.new(taxon_first_found_with_pattern.size)
		taxon_first_found_with_pattern.keys.sort.each_with_index do |taxon, ind|
			name = convert_string_to_chopped_fasta_header( taxon )
			struct = taxon_first_found_with_pattern[taxon]

			output[ind] = [name,struct].join("")
		end

		reduced_output = Sequence.remove_common_gaps(output,
			{delete_all_common_gaps: false, # remove all but one common gaps (of a series)
			ensure_common_gap_between_consecutive_non_gaps: @is_separate_introns_in_textbased_output, 
			gap_symbol: @exon_placeholder, # use placeholder as gap-symbol
			start_col: self.class.max_length_gene_name} 
		)

		return reduced_output.join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_taxonomy", exp
    	throw(:error)	
	end

	# this method is for genepainter webserver only
	# it generates a list of intron numbers that occur first in each taxon
	def export_as_taxonomy_list_of_intron_positions_per_taxon_only
		taxon_first_found_with_intron_numbers = {}
		
		sorted_intron_positions = @stats_per_intron_pos.keys.sort
		sorted_intron_positions.each_with_index do |intronpos, intron_number|
			introninfo = @stats_per_intron_pos[intronpos]
			if taxon = introninfo[:taxon_first_found] then 
				if ! taxon_first_found_with_intron_numbers[taxon] then 
					taxon_first_found_with_intron_numbers[taxon] = []
				end
				taxon_first_found_with_intron_numbers[taxon].push intron_number
			end
		end
		throw :no_taxonomy if taxon_first_found_with_intron_numbers.empty?
		
		output = Array.new(taxon_first_found_with_intron_numbers.size) 
		taxon_first_found_with_intron_numbers.keys.each_with_index do |taxon, ind|
			intron_numbers = taxon_first_found_with_intron_numbers[taxon]
			output[ind] = "#{taxon}:#{intron_numbers.join(",")}"
		end
		return output.join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_taxonomy", exp
    	throw(:error)			
	end

	# information about every intron position
	# - number of genes with this intron
	# (- first taxon this intron is assigned to and first descendants of that) 
	def export_as_statistics(exon_placeholder_output, intron_placeholder_output)
		# output consists of 
		# - reduced_genestructures and legend (to number the introns according to their position)
		# - info about each intron position

		# first part of output
		genestructures_with_legend = exon_intron_pattern_with_merged_pattern_spaces_and_intron_numbers(
			@aligned_genestructures, exon_placeholder_output, intron_placeholder_output)

		# second part of output
		sorted_intron_positions = @stats_per_intron_pos.keys.sort
		info_per_intronpos = Array.new( @stats_per_intron_pos.size + 1 ) # +1 for header
		info_per_intronpos[0] = "Intron\t# Introns"
		# does any intron have a last common ancestor?
		is_add_tax_info = @stats_per_intron_pos.collect{|k, v| v[:taxon_first_found]}.compact.any?
		if is_add_tax_info then 
			info_per_intronpos[0] += "\tLast common ancestor\tFirst unique ancestor (# Introns)"
		end
		sorted_intron_positions.each_with_index do |intronpos, intron_number|
			index_output_array = intron_number + 1

			human_readable_intron_number = Helper.ruby2human_counting(intron_number)
			n_introns = @stats_per_intron_pos[intronpos][:n_introns]

			info_per_intronpos[index_output_array] = human_readable_intron_number.to_s + "\t" + n_introns.to_s

			if taxon_first_found = @stats_per_intron_pos[intronpos][:taxon_first_found] then 
				first_descendents_with_occurences = get_first_uniq_ancestors_with_frequencies_by_genes(taxon_first_found, 
					@stats_per_intron_pos[intronpos][:genes])
				info_per_intronpos[index_output_array] += "\t" + taxon_first_found.capitalize + "\t" + first_descendents_with_occurences
			elsif is_add_tax_info
				# no taxon_first_found for this intron, but for at least one intron position
				info_per_intronpos[index_output_array] += "\tn.d.\tn.d"
			end
		end
		return [genestructures_with_legend, info_per_intronpos].join("\n")

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "export_as_statistics", exp
    	throw(:error)	
	end

	def get_first_uniq_ancestors_with_frequencies_by_genes(taxon_first_found, genes)

		descendants_with_occurence = {}
		genes.each do |gene_name|

			# just interested in taxonomy, not in exon-intron pattern: its ok to have a gene without genestructure
			gene = get_gene_obj_by_name(gene_name) 
			ind_taxon_first_found = gene.taxonomic_lineage.index(taxon_first_found)

			if ind_taxon_first_found && gene.taxonomic_lineage[ind_taxon_first_found+1] then 
				descendant = gene.taxonomic_lineage[ind_taxon_first_found+1] # lineage starts with root
				descendant = descendant.capitalize
			else
				descendant = "n.d."
			end
			if ! descendants_with_occurence[descendant] then 
				# initialize this descendant
				descendants_with_occurence[descendant] = 0
			end
			descendants_with_occurence[descendant] += 1

		end

		return prettify_first_uniq_ancestors_with_frequencies(descendants_with_occurence)
	end
	def prettify_first_uniq_ancestors_with_frequencies(ancestors_with_frequencies)
		data = []
		ancestors_with_frequencies.keys.sort.each do |ancestor|
			freq = ancestors_with_frequencies[ancestor]
			data.push( "#{ancestor} (#{freq})" )
		end
		return data.join(", ")
	end

	def exon_intron_pattern_with_merged_pattern_spaces_and_intron_numbers(not_reduced_pattern, exon_placeholder_output, intron_placeholder_output)

		pattern_names = @names_aligned_genestructures

		# add merged pattern to pattern
		# this is necessary in case all genes were analysed, but not all genes are selected for output 
		# -> intron positions that are not in patterns exist
		if @ind_merged_pattern then 
			merged_struct = @additional_structures[@ind_merged_pattern]
			merged_struct_size = merged_struct.size
		else
			# generate merged struct new
			merged_struct = generate_merged_profile(@stats_per_intron_pos, @aligned_genestructures[0].size)
			merged_struct_size = merged_struct.size
		end
		pattern_struct_size = not_reduced_pattern[0].size
		if pattern_struct_size == merged_struct_size then 
			not_reduced_pattern.push(merged_struct)
			pattern_names = [pattern_names, self.class.merged_structure_name].flatten # like this, the orginal array @names_aligned_.. is not modified
		else
			# not possible to add merged pattern
		end

		# reduce pattern (together with the merged pattern)
		reduced_pattern, reduced_merged_pattern = reduce_exon_intron_pattern( {input: not_reduced_pattern} )
		pattern = [reduced_pattern, reduced_merged_pattern].flatten

		output = Array.new(pattern.size + 1) # +1 for legend line
		n_cols = pattern.first.size 

		legend = get_empty_pattern(n_cols).split("")
		n_blanks = @stats_per_intron_pos.size.to_s.size + 1 # number of blanks needed to display intron number: the number of positions + additional blank

		# collect intronpositions in reduced structures
		intronpos_in_reduced_patterns = []
		(n_cols-1).downto(0) do |col|
			this_col = pattern.collect{ |struct| struct[col] }
			if this_col.uniq == [@exon_placeholder] then 
				# column contains only gaps
			else
				# columns contains intron
				intronpos_in_reduced_patterns.push(col)
			end
		end

		# insert gaps after each intron position in reduced structures
		pattern.each_with_index do |struct, ind|

			name = convert_string_to_chopped_fasta_header( pattern_names[ind] )
			struct = replace_exon_intron_placeholders_in_structure(struct, exon_placeholder_output, intron_placeholder_output)

			struct_with_blanks = struct.split("")

			intronpos_in_reduced_patterns.each do |pos|
				struct_with_blanks[pos] = struct_with_blanks[pos].ljust(n_blanks)
			end

			output[ind] = [name,struct_with_blanks].join("")
		end

		# generate legend with intron number
		intron_num = intronpos_in_reduced_patterns.size
		intronpos_in_reduced_patterns.each_with_index do |pos|
			legend[pos] = intron_num.to_s.ljust(n_blanks)
			intron_num -= 1
		end

		# add legend to output
		legend_name = convert_string_to_chopped_fasta_header("Intron")
		legend = legend.join("").gsub("-"," ")
		output[-1] = [legend_name,legend].join("")

		return output
	end

	# # find all introns with same phase and position
	# # depreciated, as the pure information if an intron is conserved or not is not as usefull as originally thought
	# def detect_conserved_introns(genes)
	# 	# compare introns of every gene with introns of every other genes to find duplicates
	# 	genes.combination(2) do |gene1, gene2|
	# 		common_introns = gene1.common_introns_of_this_and_other_gene(gene2)
	# 		common_introns.each do |intron|
	# 			intron.is_conserved = true
	# 		end
	# 	end

	# 	# return genes with is_conservation property set
	# 	return genes
	# end
end
