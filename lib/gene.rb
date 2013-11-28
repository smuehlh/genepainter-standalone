# a gene consists of exons and introns
class Gene
	attr_accessor :aligned_seq, :exons, :introns
	attr_reader :name

	def initialize(name)
		@name = name
		@aligned_seq = nil
		@exons = [] # exon objects in their correct order
		@introns = [] # intron objects in their correct order
	end

	def get_all_introns_with_phase
		all_pos_with_phase = @introns.collect do |intron|
			intron.get_alignmentpos_and_phase(@aligned_seq)
		end
		return all_pos_with_phase
	end

	# a conserved intron is at same position and phase as an intron in another gene
	def get_all_conserved_introns
		@introns.select do |intron|
			intron.is_conserved
		end
	end

	def get_intron_by_alignmentpos_and_phase(alignment_pos_phase)
		return @introns.find{ |i| i.get_alignmentpos_and_phase(@aligned_seq) == alignment_pos_phase }
	end

	# a sequence of same length as "@aligned_seq" consisting of gaps and intron phases only
	# exon_placeholder will be used to display exon, default: "-"
	# intron_placeholder will be used to display intron, default: intron phase
	def plot_intron_phases_onto_aligned_seq
		exon_placeholder = "-"
		intron_pos_in_alignment = Array.new(@aligned_seq.size, exon_placeholder)
		@introns.each do |intron|
			pos = intron.get_pos_in_alignment(@aligned_seq)
			intron_pos_in_alignment[pos] = intron.phase
		end

		return intron_pos_in_alignment.join("")
	end

end