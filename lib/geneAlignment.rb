# aligns gene objecs
class GeneAlignment

	# use a class level variable to make it accessible from nested class
	@@exon_placeholder = "-"
	@@intron_placeholder = nil # defaults to the intron phase
	Max_length_gene_name = 20 # also inherited, but a constant
	Struct_str = "_structure"
	# getter/setter for class variables (no instance variable, so self. is neccessary)
	def self.exon_intron_placeholder=(arr)
		@@exon_placeholder = arr[0] || "-" # default value
		@@intron_placeholder = arr[1] || nil # defaults to intron phase
	end
	def self.exon_placeholder
		@@exon_placeholder
	end
	def self.intron_placeholder
		@@intron_placeholder
	end

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
		# prepare data
		genealignment2svg_obj = GeneAlignment2svg.new(@aligned_genes, options)

		# draw genes
		output = genealignment2svg_obj.create_svg
		puts output
		
		return output
	end

	def export_as_pdb(options, consensus)
# get ref seq - check
# get seq from pdb - check
# check for special chars in either - check
# align pdb seq with referece - check
# check if alignment score is ok - check
# get introns to plot (consensus and/or of ref seq) - > use class Statistics
# plot introns onto pdb: write pymol files
		# need this to get ref seq and the introns
		alignment_genestruct = export_as_alignment_with_introns

		alignment_to_pdb_obj = GeneAlignment2pdb.new( alignment_genestruct, 
			options[:pdb_reference_protein], options[:force_alignment], options[:pdb_ref_prot_struct_only], consensus,
			options[:path_to_pdb], options[:pdb_chain], options[:pdb_penalize_endgaps]
			)
		score = alignment_to_pdb_obj.align_pdb_seq_with_ref_seq # return value 'score' is not neccessary, just to make apparent what happends
		if alignment_to_pdb_obj.is_alignment_good_enough(score) then 
			# map gene structures and write files
			alignment_to_pdb_obj.map_genestruct_onto_pdb 
		else
			Helper.warn "Cannot map gene structure onto protein structure: Alignment score is too low."
		end

		debugger
		# create obj
		# (should do alignment in init)
		# write files

		# Todo:
		# do i really need the output? maybe just a few lines to the logfile are fine ...

	end
	class Statistics
		attr_reader :n_introns_per_position

		def initialize(names_and_patterns, is_alignment)
# TODO binaere pattern, da geht addition & ist super schnell
# output merged/conserved as fasta also
			@gene_patterns = extract_exon_intron_patterns(names_and_patterns, is_alignment)
			@intron_regex = intron_placeholder_to_regexp # string or regex
			@n_introns_per_position = count_number_of_introns_per_pos_in_patterns
		end

		def intron_placeholder
			GeneAlignment.intron_placeholder
		end
		def exon_placeholder
			GeneAlignment.exon_placeholder
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