class GeneAlignment
	# aligns gene objecs
	attr_reader :aligned_genes

	def initialize(genes)
		@aligned_genes = detect_conserved_introns(genes)
	end

	# the core of this project:
	# find all introns with same phase and position
	def detect_conserved_introns(genes)

		# collect all genes (well, only the index of this gene in "genes"-Array) in which a introns with its phase occurs
		@all_conserved_introns_with_phase = {} # value: gene; key: intron+phase
		genes.each_with_index do |gene, ind|
			all_introns_with_phase = gene.get_all_introns_with_phase
			all_introns_with_phase.each do |pos_with_phase| 
				( @all_conserved_introns_with_phase[pos_with_phase] ||= [] ) << ind
			end
		end

		@all_conserved_introns_with_phase.each do |pos_phase, gene_indices|
			if gene_indices.size > 1 then
				gene_indices.each do |ind|
					requested_intron = genes[ind].get_intron_by_alignmentpos_and_phase(pos_phase)
					requested_intron.is_conserved = true
				end
			end
		end

		# return genes with is_conservation property set
		return genes
	end

	def export_as_alignment_with_introns
		output = []
		@aligned_genes.each do |gene|
			output << Sequence.convert_strings_to_fasta(gene.name, gene.aligned_seq)
			output << Sequence.convert_strings_to_fasta(gene.name + "_structure", 
				gene.plot_intron_phases_onto_aligned_seq)
		end

		return output.join("\n")
	end

	def export_as_plain_txt(exon_placeholder="-", intron_placeholder=nil)
	# TODO
	# can be faster by constructing output in look through @aligned_genes
	# than: add offset to Sequence.remove_common_gaps_from_gene_structurs(output, offset=20)
		exon_intron_patterns = []
		genenames = []
		output = [] # genename and pattern merged into one string for every genename
		@aligned_genes.each do |gene|

			exon_intron_patterns << gene.plot_intron_phases_onto_aligned_seq
			genenames << [">", gene.name.ljust(20)].join("")
		end

		reduced_exon_intron_patterns = Sequence.remove_common_gaps_from_gene_structurs(exon_intron_patterns)
		
		# interleave exon_intron_patterns with genenames
		# convert exon_intron_pattern to desired placeholders
		genenames.each_with_index do |gene, ind|
			pattern = reduced_exon_intron_patterns[ind]
			if exon_placeholder != "-" then
				pattern = pattern.gsub("-", exon_placeholder)
			end
			if intron_placeholder then
				pattern = pattern.gsub(/[0|1|2|?]/, intron_placeholder)
			end
			output << [gene,pattern].join("")
		end

		return output.join("\n")
	end

	def export_as_binary_alignment
		output = []
		exon_intron_patterns = []
		genenames = []

		@aligned_genes.each do |gene|
			exon_intron_patterns << gene.plot_intron_phases_onto_aligned_seq
			genenames << gene.name
		end

		reduced_exon_intron_patterns = Sequence.remove_common_gaps_from_gene_structurs(exon_intron_patterns)

		genenames.each_with_index do |gene, ind|
			pattern = reduced_exon_intron_patterns[ind]
			pattern = pattern.gsub(/[0|1|2|?]/, '1')
			pattern = pattern.gsub('-', '0')
			output << Sequence.convert_strings_to_fasta(gene, pattern)
		end
		
		return output.join("\n")
	end
# might need a subclass for SVG output
# svg: make gaps in exons distingishable from sequence in exons!
end