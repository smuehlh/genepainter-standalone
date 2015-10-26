class Intron
	attr_reader :pos_last_nt_in_dna_seq_before_intron, :phase, :n_nucleotides,

	# properties which cannot be set at initialization
		# pos mapped onto aligned protein sequence, but still in nucleoides instead of amino acids
		:pos_last_aa_in_aligned_protein_before_intron 
	# attr_accessor :is_conserved # true if any other intron is at same position in same phase

	def initialize(pos, length, phase)
		@pos_last_nt_in_dna_seq_before_intron = pos
		@n_nucleotides = length
		@phase = validate_intron_phase(phase) # one character string: ["0"|"1"|"2"|"?"]
	end

	# make two introns comparable
	# same if they have same position in alignment and same phase
	def eql?(other_intron)
		@pos_last_aa_in_aligned_protein_before_intron == other_intron.pos_last_aa_in_aligned_protein_before_intron &&
			@phase == other_intron.phase
	end
	def hash
		@pos_last_aa_in_aligned_protein_before_intron.hash ^ @phase.hash
	end
	# _end_ make two introns comparable

	# make two introns sortable
	def <=>(other_intron)
		[@pos_last_aa_in_aligned_protein_before_intron,@phase] <=> [
			other_intron.pos_last_aa_in_aligned_protein_before_intron, other_intron.phase] 
	end
	# _end_ make to introns sortable

	def validate_intron_phase(phase)
		phase = phase.to_s
		valid_phases = ["0", "1", "2"] # these phases are resonable
		if ! valid_phases.include?(phase) then
			phase = "?"
		end
		return phase
	end

	def self.is_phase_valid_but_non_descriptive(phase)
		return phase == "?"
	end

	def set_variables_describing_intron_in_aligned_seq(aligned_seq, is_pos_before_gap=true)
		# start and end positions are in nucleotides, aligned sequece consists of amion acids
		# 3 nucleotides code for 1 amino acid
		
		# introns of phase 0 should go underneath the first amino acid of the next exon
		# introns of phase 1/2 should go underneath the splitted amino acid
		protein_pos = convert_dna_pos_into_protein_pos(@pos_last_nt_in_dna_seq_before_intron)

		@pos_last_aa_in_aligned_protein_before_intron = 
			Sequence.sequence_pos2alignment_pos(protein_pos, aligned_seq, is_pos_before_gap) 
	end

	def convert_dna_pos_into_protein_pos(pos)
		# convert position in dna sequence into position in amino acid sequence
		# floor is important to "fix" introns of phase 2

		return (pos / 3.0).floor
	end

	def get_alignmentpos_and_phase
		return [@pos_last_aa_in_aligned_protein_before_intron, @phase]
	end

	def get_alignmentpos_merged_with_phase
		return Intron.merge_position_and_phase( *get_alignmentpos_and_phase )
	end


	# fit gene into range
	# input: n_del_nt: number of nucleotides deleted (in comparision to position of this intron)
	# input: aligned_seq: alignded sequence within range
	def create_copy_with_shifted_positions(n_del_nt, aligned_seq)
		copy = Intron.new(
			@pos_last_nt_in_dna_seq_before_intron - n_del_nt,
			@n_nucleotides,
			@phase
		)
		copy.set_variables_describing_intron_in_aligned_seq(aligned_seq)
		return copy
	end

	def self.merge_position_and_phase(pos, phase)
		pos.to_i + phase.to_i * 0.1
 	end
 	def self.are_introns_within_range(pos_phase1, pos_phase2, range)
 		diff_pos1_pos2 = 
	 		self.convert_position_and_phase_into_nt_position(pos_phase1) - self.convert_position_and_phase_into_nt_position(pos_phase2)
	 	return diff_pos1_pos2.abs <= range	
 	end

	def self.convert_position_and_phase_into_nt_position(pos_phase)
 		arr = pos_phase.divmod(1) # split float into integer and remainder
 		return arr.first * 3 + arr.last * 10 # *3 to convert into amino acid count, *10 to convert into nucleotide count
 	end

end
