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
	def initialize(genes, consensus_val, is_merged_pattern, is_statistics, taxonomy_options)
		# @aligned_genes = detect_conserved_introns(genes)
		@genes = genes

		@exon_placeholder = "-"
		@intron_placeholder = nil # defaults to intron phase

		# overwritten by method convert_to_exon_intron_pattern
		@ind_consensus_pattern = nil
		@ind_merged_pattern = nil
		@ind_tax_pattern = nil
		@stats_per_intron_pos = {}
		
		# convert_to_exon_intron_pattern overwrites also @ind_consensus_pattern, @ind_merged_pattern, @ind_tax_pattern if neccessary
		@aligned_genestructures = convert_to_exon_intron_pattern(consensus_val, is_merged_pattern, is_statistics, taxonomy_options) # are of same order as @genes !!!
		@reduced_aligned_genestructures = reduce_exon_intron_pattern # reduce gene structures to "needed" parts 

		@n_structures = @aligned_genestructures.size
	end

	# align genestructures by plotting them onto the multiple sequence alignment
	# exon representation: "-", intron representation: phase
	# output
	# Array of strings: exon intron pattern plotted onto the aligned sequence
	def convert_to_exon_intron_pattern(consensus_val, is_merged_pattern, is_statistics, taxonomy_options)
		# collect the exon_intron_patterns of each gene, also for merged/consensus

		patterns = Array.new(@genes.size)

		# number of introns, phase, and genes it occurs in for each intron position
		intron_pos_with_annotation = Hash.new { |h,k| h[k] = { n_introns: 0, phase: "?", genes: [] } }

		# parse taxonomy input
		genes_within_selected_taxa = taxonomy_options[:genes_within_taxa] || []
		is_intron_exclusive_for_taxa = taxonomy_options[:is_exclusive]

		@genes.each_with_index do |gene, ind|
			patterns[ind] = gene.plot_intron_phases_onto_aligned_seq(@exon_placeholder, @intron_placeholder)

			update_annotation_per_intron_position(gene, intron_pos_with_annotation)
		end

		if consensus_val then 
			pattern_length = patterns[0].size
			min_n_introns_per_pos = consensus_val * @genes.size
			patterns << generate_consensus_profile( intron_pos_with_annotation, min_n_introns_per_pos , pattern_length )
			@ind_consensus_pattern = patterns.size - 1 # indices start with 0
			self.class.consensus_val = consensus_val # set class variable needed for export functions
		end
		if is_merged_pattern then
			pattern_length = patterns[0].size
			patterns << generate_merged_profile( intron_pos_with_annotation, pattern_length )
			@ind_merged_pattern = patterns.size - 1 # indices start with 0
		end
		if genes_within_selected_taxa.any? then
			pattern_length = patterns[0].size
			patterns << generate_taxonomy_pattern( intron_pos_with_annotation, genes_within_selected_taxa, is_intron_exclusive_for_taxa, 
				pattern_length )
			@ind_tax_pattern = patterns.size - 1 # indices start with 0
		end
		if is_statistics then 
			# stats also include phase information, which are currently not used
			@stats_per_intron_pos = intron_pos_with_annotation
		end

		return patterns
	end

	# methods for "statistics per intron position and statistics: merged/consensus/taxonomy profile"

	# update counts, phase and occurence-list per intron position
	def update_annotation_per_intron_position(gene, intron_pos_with_annotation)
		gene.get_all_introns_with_phase.each do |pos, phase|

			intron_pos_with_annotation[pos][:n_introns] += 1
			intron_pos_with_annotation[pos][:phase] = phase
			intron_pos_with_annotation[pos][:genes].push gene.name
		end
	end
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
	def get_empty_pattern(len)
		return @exon_placeholder * len
	end
	def generate_taxonomy_pattern( intron_pos_with_annotation, genes_within_selected_taxa, is_intron_exclusive_for_taxa, pattern_length )
		tax_pattern = get_empty_pattern(pattern_length)

		intron_pos_with_annotation.each do |intronpos, info|
			if ( is_intron_exclusive_for_taxa && (genes_within_selected_taxa == info[:genes]) ) ||
				( ! is_intron_exclusive_for_taxa && (genes_within_selected_taxa & info[:genes]).any? ) then
				tax_pattern[intronpos] = info[:phase]
			end
		end
		return tax_pattern
	end
	# end methods for statistics

	def reduce_exon_intron_pattern
		# important: duplicate the pattern ! 
		patterns = @aligned_genestructures.map {|ele| ele.dup}
		Sequence.remove_common_gaps(patterns,
			{is_alignment: false, # input is not fasta-formatted alignment
			gap_symbol: @exon_placeholder} # use placeholder as gap-symbol
			)
	end

	def convert_string_to_chopped_fasta_header(str)
		# add ">" at beginning of string
		# left justify string to certain length
		(">" << str).ljust(self.class.max_length_gene_name)
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

		@aligned_genestructures.each_with_index do |struct, ind|

			index_output_struct = ind * 2 + 1
			if ind == @ind_merged_pattern then 
				name_struct = self.class.merged_structure_name
			elsif ind == @ind_consensus_pattern
				name_struct = self.class.consensus_structure_name
			elsif ind == @ind_tax_pattern
				name_struct = self.class.taxonomy_structure_name
			else
				# its the structure of a @aligned_gene
				index_output_seq = ind * 2
				gene = @genes[ind]
				output[index_output_seq] = Sequence.convert_strings_to_fasta(gene.name, gene.aligned_seq)
				name_struct = gene.name + self.class.suffix_structure_in_alignment
			end 

			output[index_output_struct] = Sequence.convert_strings_to_fasta(name_struct, struct)
		end

		return output.join("\n")
	end

	def export_as_binary_alignment
		# output contains binary gene structure, in fasta format
		output = Array.new(@n_structures)

		exon_placeholder_output = "0"
		intron_placeholder_output = "1"

		@reduced_aligned_genestructures.each_with_index do |struct, ind|
			if ind == @ind_merged_pattern then 
				name = self.class.merged_structure_name
			elsif ind == @ind_consensus_pattern 
				name = self.class.consensus_structure_name
			elsif ind == @ind_tax_pattern
				name = self.class.taxonomy_structure_name
			else
				# its a gene
				name = @genes[ind].name
			end
			struct = replace_exon_intron_placeholders_in_structure(struct, exon_placeholder_output, intron_placeholder_output)

			output[ind] = Sequence.convert_strings_to_fasta( name, struct )		
		end

		return output.join("\n")
	end

	def export_as_plain_txt(exon_placeholder_output, intron_placeholder_output)
	
		output = Array.new(@n_structures)

		@reduced_aligned_genestructures.each_with_index do |struct, ind|
			if ind == @ind_merged_pattern then 
				name = self.class.merged_structure_name
			elsif ind == @ind_consensus_pattern
				name = self.class.consensus_structure_name
			elsif ind == @ind_tax_pattern
				name = self.class.taxonomy_structure_name
			else
				# its a gene
				name = @genes[ind].name
			end

			name = convert_string_to_chopped_fasta_header ( name )
			struct = replace_exon_intron_placeholders_in_structure( struct, exon_placeholder_output, intron_placeholder_output )

			output[ind] = [name, struct].join("")
					
		end

		return output.join("\n")
	end

	def export_as_svg(options)
		# prepare data
		genealignment2svg_obj = GeneAlignment2svg.new(@genes, options)

		# draw genes
		output = genealignment2svg_obj.create_svg
		
		return output
	end

	def export_as_pdb(options)

		output = []
		Helper.log "\nMap gene structures onto PDB #{options[:path_to_pdb]}"
		GeneAlignment2pdb.log_howtouse_pythonscripts

		# prepare data
		ref_seq = "" # the aligned reference sequence
		ref_struct = "" # the genestructure plotted onto the ref_seq

		# default: use first seq in alignment as reference
		# if someother gene is specified via command line, overwrite sequence index
		ind_of_ref_seq = 0
		if options[:pdb_reference_protein] then 
			# use specified seq
			name = options[:pdb_reference_protein]
			if name.starts_with?(">") then
				name = name[1..-1]
			end
			ind_of_ref_seq = @genes.collect{ |g| g.name }.index(name)
			if ! ind_of_ref_seq then 
				Helper.log "Cannot match gene #{options[:pdb_reference_protein]}."
				Helper.log "Use the first gene in alignment instead."
				ind_of_ref_seq = 0
			end
		end

		ref_gene = @genes[ind_of_ref_seq]
		ref_seq = ref_gene.aligned_seq

		# get gene structure
		if options[:pdb_ref_prot_struct_only]
			# structure of the reference sequence only
			ind_ref_struct = ind_of_ref_seq
		else
			# need complete alignment with intron pos and phases
			# to get merged or consensus sequence
			if @ind_merged_pattern then 
				# use the merged pattern
				ind_ref_struct = @ind_merged_pattern
			end
			if @ind_consensus_pattern then 
				ind_ref_struct = @ind_consensus_pattern
			end
			if @ind_tax_pattern then 
				ind_ref_struct = @ind_tax_pattern
			end
		end
		ref_struct = @aligned_genestructures[ind_ref_struct]

		genealignment2pdb_obj = GeneAlignment2pdb.new( ref_seq, ref_struct, options )

		# plot gene structures onto pdb
		output1, output2 = genealignment2pdb_obj.map_genestructure_onto_pdb

		genealignment2pdb_obj.log_alignment(ref_gene.name, options[:path_to_pdb])

		return output1, output2
	end

	def export_as_taxonomy(bla)
		# compare consensus pattern of every node with the consensus pattern of every node before (direction to root)

		# wie zusammenspiel mit anderen tax -optionen???

	end

	def export_as_statistics(ancestors_of_each_gene)

		# output contains of
		# 1) exon-intron patterns which labels each intron position with an number
		# 2) information about each intron position: number of introns , last common ancestor, ancestor following the LCA

		# collect data for output
		# => for each pattern:  pattern with blanks & legend-like lines
		genestructures_with_legend = exon_intron_pattern_with_spaces_and_legend_for_intron_numbers

		# => for each intron position: lca, list with first ancestors after lca with n_genes having this
		infos_per_intronpos = Array.new( @stats_per_intron_pos.size  + 1)
		infos_per_intronpos[0] = "Intron number\t # introns\t last comman ancestor; first ancestors after the last common ancestor (if applicable)"
		@stats_per_intron_pos.keys.sort.each_with_index do |intron_pos, ind|
			# ind is number of intron, not number of introns occuring at this position

			ind_info_arr = ind + 1
			one_based_ind = Helper.ruby2human_counting(ind) # by change, this is the same as ind_info_arr ... but its cleaner with 2 vars

			n_introns = @stats_per_intron_pos[intron_pos][:n_introns]

			genes = @stats_per_intron_pos[intron_pos][:genes]
			if ancestors_of_each_gene then 
				ancestors_of_these_genes = genes.collect{ |g| ancestors_of_each_gene[g] }
				ancestors_of_these_genes.delete_if{|sub_arr| sub_arr == [] }
			else
				ancestors_of_these_genes = []
			end

			if ancestors_of_these_genes.any? then 
				last_common = Taxonomy.last_common_ancestor_of_genes(ancestors_of_these_genes)
				first_uniq_list = Taxonomy.first_uniq_ancestors_of_genes(ancestors_of_these_genes)
				first_uniq_freq = Helper.word_frequency(first_uniq_list)
				first_uniq_pretty_arr = []
				# sort frequencies in reverse order
				first_uniq_freq.sort_by {|_k,v| v}.reverse.each do |word, freq|
					first_uniq_pretty_arr << "#{word} (#{freq})"
				end

				# save infos for this intron position
				infos_per_intronpos[ind_info_arr] = "#{one_based_ind}\t#{n_introns}\t#{last_common}; #{first_uniq_pretty_arr.join(", ")}"
			else
				infos_per_intronpos[ind_info_arr] = "#{one_based_ind}\t#{n_introns}\tNo information about taxonomy."
			end

		end

		return [ genestructures_with_legend, infos_per_intronpos ].join("\n")

	end

	def exon_intron_pattern_with_spaces_and_legend_for_intron_numbers

		output = Array.new(@genes.size + 1) # +1 for legend-like line
		sorted_intronpositions = @stats_per_intron_pos.keys.sort
		n_blanks = sorted_intronpositions.size.to_s.size + 1 # number of blanks needed to display intron number: the number of positions + additional blank
		legend = get_empty_pattern(@aligned_genestructures.first.size).split("")

		# iterate over genes instead of aligned_genestructures to avoid merged/consensus/taxonomy pattern
		@genes.each_with_index do |gene, ind|

			name = convert_string_to_chopped_fasta_header( gene.name )
			struct = @aligned_genestructures[ind]

			# insert blanks after each intron position
			struct = replace_exon_intron_placeholders_in_structure( struct, "-", "|" )
			struct = struct.split("")

			sorted_intronpositions.each_with_index do |pos, num|
				# num is number of intron, not number of introns occuring at this position
				one_based_count = Helper.ruby2human_counting(num) 

				struct[pos] = struct[pos].ljust(n_blanks)
				legend[pos] = one_based_count.to_s.ljust(n_blanks)

			end

			output[ind] = [name, struct].join("")
		end

		output[-1] = [ convert_string_to_chopped_fasta_header( "Intron number" ), legend ].join("")
		# remove common gaps from structures
		output = Sequence.remove_common_gaps(output, {is_alignment: false})

		# replace "-" with " " to prettify the legend-string
		output[-1] = output[-1].gsub("-", " ")

		return output.join("\n")
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