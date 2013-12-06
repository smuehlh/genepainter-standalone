# a gene consists of exons and introns
class Gene
	attr_accessor :aligned_seq, :exons, :introns, :alignment_range
	attr_reader :name

	def initialize(name)
		@name = name
		@aligned_seq = nil
		@exons = [] # exon objects in their correct order
		@introns = [] # intron objects in their correct order
		@alignment_range = [] # array of arrays; restrict mapping of gene stucture onto alignment to these parts of the alignment
	end

	# set instance variable @aligned_seq and also intron position in alignment
	def add_aligned_seq(aligned_seq)
		@aligned_seq = aligned_seq
		@introns.each do |intron|
			intron.pos_last_aa_in_aligned_seq_before_intron = intron.get_pos_in_alignment(aligned_seq)
		end
	end

	def add_alignment_range(ranges)
		# make sure that the very last position of ranges does not exceed aligned sequence size
		if ranges[-1] && ! @aligned_seq[ranges[-1][1]] then
			Helper.abort "Invalid range: Range exceeds alignment size"
		end
		@alignment_range = ranges
	end

	# returns common introns of both genes, thus, returned introns are not uniq
	def common_introns_of_this_and_other_gene(other_gene)
		(@introns & other_gene.introns).concat(other_gene.introns & @introns)
	end

	def get_aligned_seq_within_range
		reduce_str_to_range(@aligned_seq)
	end

	def get_all_introns_with_phase
		all_pos_with_phase = @introns.collect do |intron|
			intron.get_alignmentpos_and_phase
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
	# exon_representation will be used to display exon, default: "-"
	# intron_representation will be used to display intron, default: intron phase
	def plot_intron_phases_onto_aligned_seq(exon_representation="-", intron_representation=nil)

		intron_pos_in_alignment = Array.new(@aligned_seq.size, exon_representation)
		@introns.each do |intron|
			pos = intron.pos_last_aa_in_aligned_seq_before_intron
			intron_pos_in_alignment[pos] = intron_representation || intron.phase
		end

		reduced_pattern = reduce_str_to_range( intron_pos_in_alignment.join("") )
		return reduced_pattern
	end

	# reduce the aligned sequence or the exon_intron_pattern to @alignment_range
	def reduce_str_to_range(str)
		reduced_str = nil # will be converted to empty string if there are any ranges
		@alignment_range.each do |start_at, stop_at|
			# collect all the needed parts of str
			reduced_str = reduced_str.to_s << str.slice(start_at..stop_at)
		end
		return reduced_str || str # if no ranges, str will be returned
	end

end