class Exon
	attr_reader :start_pos_in_dna_seq, :end_pos_in_dna_seq,

		# these attributes cannot be set at initialization, as they require the aligned sequence
		:start_pos_in_aligned_protein, :end_pos_in_aligned_protein

	def initialize(start, stop)
		@start_pos_in_dna_seq = start
		@end_pos_in_dna_seq = stop
	end

	def length_in_nt
		@end_pos_in_dna_seq - @start_pos_in_dna_seq
	end

	def length_in_alignment
		@end_pos_in_aligned_protein - @start_pos_in_aligned_protein
	end

	def set_variables_describing_exon_in_aligned_seq(aligned_seq)
		@start_pos_in_aligned_protein = 
			Sequence.sequence_pos2alignment_pos(convert_dna_pos_into_protein_pos(@start_pos_in_dna_seq), aligned_seq)
		@end_pos_in_aligned_protein = 
			# the last position must be excluded, it belongs actually to the next exon
			Sequence.sequence_pos2alignment_pos(convert_dna_pos_into_protein_pos(@end_pos_in_dna_seq) - 1, aligned_seq)
	end

	def convert_dna_pos_into_protein_pos(pos_in_dna)
		# start and end positions are in nucleotides, aligned sequece consists of amion acids
		# 3 nucleotides code for 1 amino acid

		# CAUTION: this method does not translate into the position in the aligned sequence (there would be gaps also!)
		(pos_in_dna/3.0).floor
	end


	# fit gene into range
	# input: n_del_nt: number of nucleotides deleted (in comparision to position of this exon)
	# input: start_pos_nt: start position in dna sequence [ignored if nil]
	# input: end_pos_nt: end position in dna sequence [ignored if nil]
	# input: aligned_seq: alignded sequence within range
	def create_copy_with_shifted_positions(n_del_nt, start_in_dna, end_in_dna, aligned_seq)
		copy = Exon.new(
			start_in_dna || @start_pos_in_dna_seq - n_del_nt, # start position, reuse start of this exon if possible
			end_in_dna || @end_pos_in_dna_seq - n_del_nt # end position, reuse end of this exon if possible
			)
		copy.set_variables_describing_exon_in_aligned_seq(aligned_seq)
		return copy
	end

end
