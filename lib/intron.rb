class Intron
	attr_reader :pos_last_aa_in_protein_seq_before_intron, :phase, :n_nucleotides

	# properties which cannot be set at initialization
	attr_accessor :is_conserved, # true if any other intron is at same position in same phase
		:pos_last_aa_in_aligned_seq_before_intron # pos mapped onto aligned protein sequence

	def initialize(pos, length, phase)
		@pos_last_aa_in_protein_seq_before_intron = pos - 1
		@n_nucleotides = length
		@phase = validate_intron_phase(phase) # one character string: ["0"|"1"|"2"|"?"]
		@is_conserved = false # "default" value is very pessimistic
		@pos_last_aa_in_aligned_seq_before_intron = nil
	end

	# make two introns comparable
	# same if they have same position in alignment and same phase
	def eql?(other_intron)
		@pos_last_aa_in_aligned_seq_before_intron == other_intron.pos_last_aa_in_aligned_seq_before_intron &&
			@phase == other_intron.phase
	end
	def hash
		@pos_last_aa_in_aligned_seq_before_intron.hash ^ @phase.hash
	end
	# _end_ make two introns comparable

	def validate_intron_phase(phase)
		phase = phase.to_s
		valid_phases = ["0", "1", "2"] # these phases are resonable
		if ! valid_phases.include?(phase) then
			phase = "?"
		end
		return phase
	end

	def is_phase_valid_but_non_descriptive(phase)
		return phase == "?"
	end

	def get_pos_in_alignment(aligned_seq)
		return Sequence::sequence_pos2alignment_pos(@pos_last_aa_in_protein_seq_before_intron, aligned_seq)
	end

	def get_alignmentpos_and_phase
		return [@pos_last_aa_in_aligned_seq_before_intron, @phase]
	end

end