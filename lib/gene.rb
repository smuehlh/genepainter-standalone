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

	def get_all_gap_pos
		Helper.find_each_index(@aligned_seq, "-")
	end

	def get_all_gaps_with_length(is_convert_to_nt_length=false)
		all_gaps = get_all_gap_pos
		if is_convert_to_nt_length then
			length_one_gap = 3
		else
			length_one_gap = 1
		end

		# initialize: the very first pos is a gap of minimum length_one_gap
		all_pos_with_length = [[all_gaps[0],length_one_gap]]
		all_gaps.each_cons(2) do |x,y|
			if y == x+length_one_gap then
				# same gap, make it a bit longer
				all_pos_with_length[-1][1] += length_one_gap
			else
				# new gap
				all_pos_with_length << [y*length_one_gap,length_one_gap]
			end
		end
		return all_pos_with_length
	end

	def get_all_exons_with_length(is_convert_to_nt_length=false)
		if is_convert_to_nt_length then
			length_one_pos = 3
		else
			length_one_pos = 1
		end		
		all_pos_with_length = @exons.collect do |exon|
			[exon.get_start_pos_in_alignment(@aligned_seq) * length_one_pos , exon.aligned_seq_length(@aligned_seq) * length_one_pos ]
		end
		return all_pos_with_length
	end

	def get_all_introns_with_length(is_convert_to_nt_length=false)
		if is_convert_to_nt_length then
			length_one_pos = 3
		else
			length_one_pos = 1
		end	
		all_pos_with_length = @introns.collect do |intron|
			[intron.pos_last_aa_in_aligned_seq_before_intron * length_one_pos, intron.n_nucleotides]
		end
		return all_pos_with_length
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
		@introns.find{ |i| i.get_alignmentpos_and_phase(@aligned_seq) == alignment_pos_phase }
	end

	def length_of_exons_in_aa
		sum = 0
		@exons.each { |exon| sum += (exon.end_pos_in_protein_seq - exon.start_pos_in_protein_seq) }
		return sum.to_f
	end

	def length_of_introns_in_nt
		sum = 0
		@introns.each { |intron| sum += intron.n_nucleotides }
		return sum.to_f
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