class Exon
	attr_reader :start_pos_in_protein_seq, :end_pos_in_protein_seq

	def initialize(start, stop)
		@start_pos_in_protein_seq = start
		@end_pos_in_protein_seq = stop

		# TODO
		# also needs: pos_sequence_shifts = [] ?
		# also needs: n_amino_acids = endpos-startpos ?
	end

	def get_start_pos_in_alignment(aligned_seq)
		return Sequence::sequence_pos2alignment_pos(@start_pos_in_protein_seq, aligned_seq)
	end

	def get_end_pos_in_alignment(aligned_seq)
		return Sequence::sequence_pos2alignment_pos(@end_pos_in_protein_seq, aligned_seq)
	end

	def get_all_gap_pos_in_alignment(aligned_seq)
		gap_symbol = "-" 
		pos = []
		i = -1 # offset in aligned_seq for search for occurance
		while i = aligned_seq.index(gap_symbol, i+1)
			pos << i
		end
		return pos
	end

end