class Exon
	attr_reader :start_pos_in_protein_seq, :end_pos_in_protein_seq

	def initialize(start, stop)
		@start_pos_in_protein_seq = start
		@end_pos_in_protein_seq = stop

		# TODO
		# also needs: pos_sequence_shifts = [] ?
	end

	def get_start_pos_in_alignment(aligned_seq)
		return Sequence.sequence_pos2alignment_pos(@start_pos_in_protein_seq, aligned_seq)
	end

	def get_end_pos_in_alignment(aligned_seq)
		return Sequence.sequence_pos2alignment_pos(@end_pos_in_protein_seq, aligned_seq)
	end

	def aligned_seq_length(aligned_seq)
		return get_end_pos_in_alignment(aligned_seq) - get_start_pos_in_alignment(aligned_seq)
	end
end