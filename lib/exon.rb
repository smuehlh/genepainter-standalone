class Exon
	attr_reader :start_pos_in_dna_seq, :end_pos_in_dna_seq,

		# these attributes cannot be set at initialization, as they require the aligned sequence
		:start_pos_in_aligned_protein, :end_pos_in_aligned_protein, 
		:length_in_alignment

	def initialize(start, stop)
		@start_pos_in_dna_seq = start
		@end_pos_in_dna_seq = stop

		# TODO
		# also needs: pos_sequence_shifts = [] ?
	end

	def set_variables_describing_exon_in_aligned_seq(aligned_seq)
		@start_pos_in_aligned_protein = 
			Sequence.sequence_pos2alignment_pos(convert_dna_pos_into_protein_pos(@start_pos_in_dna_seq), aligned_seq)
		@end_pos_in_aligned_protein = 
			Sequence.sequence_pos2alignment_pos(convert_dna_pos_into_protein_pos(@end_pos_in_dna_seq), aligned_seq)

		@length_in_alignment = @end_pos_in_aligned_protein - @start_pos_in_aligned_protein
	end

	def convert_dna_pos_into_protein_pos(pos_in_dna)
		# start and end positions are in nucleotides, aligned sequece consists of amion acids
		# 3 nucleotides code for 1 amino acid
		(pos_in_dna/3.0).round
	end

end
