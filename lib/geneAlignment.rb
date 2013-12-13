# aligns gene objecs
class GeneAlignment

	attr_reader :aligned_genes

	# class variables are shared by subclass Statistics
	@@exon_placeholder = "-"
	@@intron_placeholder = nil # defaults to the intron phase
	Max_length_gene_name = 20 # also inherited, but a constant
	Struct_str = "_structure"

	def initialize(genes)
		@aligned_genes = detect_conserved_introns(genes)
	end

	# find all introns with same phase and position
	def detect_conserved_introns(genes)
		# compare introns of every gene with introns of every other genes to find duplicates
		genes.combination(2) do |gene1, gene2|
			common_introns = gene1.common_introns_of_this_and_other_gene(gene2)
			common_introns.each do |intron|
				intron.is_conserved = true
			end
		end

		# return genes with is_conservation property set
		return genes
	end

	# align genestructures by plotting them onto the MSA and reducing the gene structures to common gaps
	# output format: array of strings, containing concatenated gene name and reduced exon intron pattern
	# exons = "-", introns = phase
	def convert_to_exon_intron_pattern
		# collect gene names and the exon_intron_patterns
		names_and_patterns = Array.new(@aligned_genes.size) # genename and exon-intron pattern
		@aligned_genes.each_with_index do |gene, ind|
			name = (">" << gene.name).ljust(Max_length_gene_name)
			exon_intron_pattern = gene.plot_intron_phases_onto_aligned_seq(@@exon_placeholder, @@intron_placeholder)
			names_and_patterns[ind] = [name, exon_intron_pattern].join("")
		end

		# remove common gaps from exon_intron_patterns
		names_and_reduced_patterns = Sequence.remove_common_gaps(names_and_patterns, 
			false,  # false: input is not a sequence alignment
			Max_length_gene_name, @@exon_placeholder)

		return names_and_reduced_patterns
	end

	def export_as_alignment_with_introns
		output = []
		@aligned_genes.each do |gene|
			output << Sequence.convert_strings_to_fasta(gene.name, gene.get_aligned_seq_within_range)
			output << Sequence.convert_strings_to_fasta(gene.name + Struct_str, 
				gene.plot_intron_phases_onto_aligned_seq)
		end

		return output
	end

	def export_as_plain_txt
		output = convert_to_exon_intron_pattern

		return output
	end

	def export_as_binary_alignment
		output = convert_to_exon_intron_pattern

		# convert into fasta-formatted string
		converted_output = output.map do |name_and_pattern|
			name, pattern = name_and_pattern[0..Max_length_gene_name-1], name_and_pattern[Max_length_gene_name..-1]
			name_and_pattern = Sequence.convert_strings_to_fasta(name, pattern)
		end
		
		return converted_output
	end

def export_as_svg(options)
	# the general idea:
	# need to find the longest intron at every position and expand the shorter
	# and have a method of Svg for every type of variable, such that the Svg method does not need to know that much about the underlying structure!!!

	genes_with_aligned_features = {}

	@aligned_genes.combination(2) do |gene1, gene2|

		## prepare new data structure for genes (if neccessary)
		if ! genes_with_aligned_features[gene1.name] then
			genes_with_aligned_features[gene1.name] = PrepareGeneForSvg.new_data_struct_for_gene(gene1)
		end
		if ! genes_with_aligned_features[gene2.name] then 
			genes_with_aligned_features[gene2.name] = PrepareGeneForSvg.new_data_struct_for_gene(gene2)
		end

		## are introns common to both genes?
		# common_introns contain introns of both genes
		common_introns = gene1.common_introns_of_this_and_other_gene(gene2)
		## no: nothing to do
		## yes: iterate through every pair of common introns
		common_introns.sort.each_slice(2) do |intron1, intron2|
			# find out which intron belongs to which gene 
			if gene1.introns.include?(intron1) then
				gene_features_belonging_to_intron1 = genes_with_aligned_features[gene1.name]
				gene_features_belonging_to_intron2 = genes_with_aligned_features[gene2.name]
				ind_intron1 = gene1.introns.index(intron1)
				ind_intron2 = gene2.introns.index(intron2)
			else
				gene_features_belonging_to_intron1 = genes_with_aligned_features[gene2.name]
				gene_features_belonging_to_intron2 = genes_with_aligned_features[gene1.name]
				ind_intron1 = gene2.introns.index(intron1)
				ind_intron2 = gene1.introns.index(intron2)
			end

			## 	compare length of g1 intron and g2 intron
			## the gene with the shorter intron needs some more offset to keep consequtive features aligned
			length_diff_intron1_intron2 = intron1.n_nucleotides - intron2.n_nucleotides
			if length_diff_intron1_intron2 < 0 then
				# intron 1 is shorter
				PrepareGeneForSvg.adding_offset(gene_features_belonging_to_intron1, length_diff_intron1_intron2, ind_intron1)
			else
				# intron 2 is shorter
				PrepareGeneForSvg.adding_offset(gene_features_belonging_to_intron2, length_diff_intron1_intron2, ind_intron2)
			end
		end

	end
debugger
puts "next: scaling aligned features to meet ratio exons & introns and the total width"
	# TODO
	# add offsets to positions
	# scale features (maybe collect some information for scaling during the processing, like the max lenght of exons and introns found in any gene)
	# draw featrures 
	# -(rewrite Svg: make one method to draw boxes and instruct it with the size and the color)
	# -every type of feature should be printed in one rush (at moment, its more like every gene is printed in one rush, which require lots of logic about feature in Svg)

	# idea for 'reduced': replace intron length by some fixed number!

	# idea for drawing: start with large box (balken im hintergrund) for intron-extension
	# plot boxes for features (exons, gaps and real introns) onto this 
end

	def export_as_svg_alter_ansatz(options)
		output = []
		n_genes = @aligned_genes.size

		## parameters for drawing
		# params = Svg.parameters
		size_per_gene = Svg.size_per_gene(n_genes) # height of each gene, including the space between genes
		# height_per_gene_feature = Svg.height_of_gene_feature(size_per_gene) # size of boxes representing exons and introns

		## collecing data
		names = []
		genes = []
		y_pos = []

		longest_aligned_seq = 0 # in nucleotides, counted per gene
		longest_intronic_seq = 0 # in nucleotides, counted per gene
		@aligned_genes.each_with_index do |gene, ind|

			this_length_introns = gene.length_of_introns_in_nt
			this_length_aligned_exons = gene.aligned_seq.size * 3

			if this_length_introns > longest_intronic_seq then
				longest_intronic_seq = this_length_introns
			end
			if this_length_aligned_exons > longest_aligned_seq then
				longest_aligned_seq = this_length_aligned_exons
			end

			genes << {
				exons: gene.get_all_exons_with_length, 
				gaps: gene.get_all_gaps_with_length,
				introns: gene.get_all_introns_with_length
			}
			names << gene.name[0..Max_length_gene_name]
			y_pos << (ind * size_per_gene)
		end
# TODO
# same color for conserved introns
# via stats: assign color based on pos! and than count introns before to know which intron gets which color
		## convert data to picture
		scaled_genes = Svg.scale_genes_to_canvas(genes, longest_aligned_seq, longest_intronic_seq)

		output << Svg::Painter.header(options[:size][:width], options[:size][:height])
		output << Svg.print_names_and_genes(names, scaled_genes, y_pos, options[:reduced]) 
		output << Svg::Painter.footer
debugger
		return output
	end



# might need a subclass for SVG output
# svg: make gaps in exons distingishable from sequence in exons!
# svg: make use of length stored in every exon and intron object
# and also (for coloring the simple output): of the property: _is_conserved!!!

	module PrepareGeneForSvg
		def self.new_data_struct_for_gene(gene)
			# each value is array of arrays (for each feature: position in alignment, offset, length)
			return {
				exons: expand_array_for_offset(gene.get_all_exons_with_length(true)), # true: convert amino acid to nucleotide count
				exon_gaps: expand_array_for_offset(gene.get_all_gaps_with_length(true)), # true: convert aa to nt count 
				introns: expand_array_for_offset(gene.get_all_introns_with_length)#,
				# intron_gaps: []
			}
		end
		def self.expand_array_for_offset(array)
			array.map! do |ele|
				ele = [ele[0],0,ele[1]]
			end	
		end
		def self.adding_offset(struct, offset, ind)
			struct.each do |key, val|
				# add offset to every value consecutive to ind
				val[ind+1..-1].map!{ |arr| arr[1] += offset }
			end
		end
	end

	class Statistics

		def initialize(names_and_patterns, is_alignment)
# TODO binaere pattern, da geht addition & ist super schnell
# output merged/conserved as fasta also
			@gene_patterns = extract_exon_intron_patterns(names_and_patterns, is_alignment)
			@intron_regex = intron_placeholder_to_regexp # string or regex
			@n_introns_per_position = count_number_of_introns_per_pos_in_patterns
		end

		def intron_placeholder
			GeneAlignment.class_variable_get(:@@intron_placeholder)
		end
		def exon_placeholder
			GeneAlignment.class_variable_get(:@@exon_placeholder)
		end

		def intron_placeholder_to_regexp
			if intron_placeholder then 
				return Regexp.new( Regexp.escape( intron_placeholder ) ) # escaping is necessary because placeholder might be "|"
			else
				# placeholder = 'nil' actually means, that intron is represented as [0|1|2|?]
				return Regexp.new( "[0|1|2|?]" )
			end
		end

		def extract_exon_intron_patterns(input, is_alignment)
			if is_alignment then
				names, seqs_and_patterns = Sequence.convert_fasta_array_back_to_arrays(input)
				if exon_placeholder == "0" && intron_placeholder == "1" then
					# gene structure in binary fasta format
					return seqs_and_patterns
				else
					# alignment with gene structures
					# return only gene structures
					return seqs_and_patterns.select.with_index { |seq, ind| seq if names[ind].include?(Struct_str) }
				end
			else
				return input.collect{|line| line[Max_length_gene_name..-1] }
			end
		end

		def count_number_of_introns_per_pos_in_patterns
			counts = Hash.new(0) # keys: position in gene_patterns, values: number of occurances
			@gene_patterns.each do |pattern|
				inds = Helper.find_each_index(pattern,@intron_regex)
				inds.each do |ind|
					counts[ind] += 1
				end
			end
			return counts
		end

		def get_consensus_exon_intron_pattern(pct_value, is_alignment=false)
			consensus_pattern = get_an_empty_pattern
			cons_name = get_name_for_new_pattern("Consensus")
			n_genes = @gene_patterns.size
			@n_introns_per_position.each do |pos, num|
				# check if the intron occurs often enough
				if num.to_f/n_genes >= pct_value then
					consensus_pattern[pos] = intron_placeholder
				end
			end
			if is_alignment then
				return Sequence.convert_strings_to_fasta(cons_name, consensus_pattern)
			else
				return (cons_name << consensus_pattern)
			end
		end

		def get_merged_exon_intron_pattern(is_alignment=false)
			merged_pattern = get_an_empty_pattern
			merged_name = get_name_for_new_pattern("Merged")
			@n_introns_per_position.keys.each do |pos|
				merged_pattern[pos] = intron_placeholder
			end
			if is_alignment then
				return Sequence.convert_strings_to_fasta(merged_name, merged_pattern)
			else
				return (merged_name << merged_pattern)
			end
		end

		def get_an_empty_pattern
			pattern_length = @gene_patterns[0].size
			new_pattern = exon_placeholder * pattern_length # 'empty' pattern, only exons
			return new_pattern
		end

		def get_name_for_new_pattern(str)
			return str.ljust(Max_length_gene_name)
		end

	end
end