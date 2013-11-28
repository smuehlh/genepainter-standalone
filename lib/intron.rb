class Intron
	attr_reader :pos_last_aa_in_protein_seq_before_intron, :phase, :n_nucleotides
	# this property will be set by class GeneAlignment
	attr_accessor :is_conserved # true if any other intron is at same position in same phase

	def initialize(pos, length, phase)
		@pos_last_aa_in_protein_seq_before_intron = pos - 1
		@n_nucleotides = length
		@phase = validate_intron_phase(phase) # one character string: ["0"|"1"|"2"|"?"]
		@is_conserved = false # "default" value is very pessimistic
	end

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

	def get_alignmentpos_and_phase(aligned_seq)
		return [get_pos_in_alignment(aligned_seq), @phase]
	end

end