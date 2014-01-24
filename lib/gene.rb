# a gene consists of exons and introns
# Important: all get_... methods should return features _within_ alignment range

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
		@exons.each do |exon|
			exon.set_variables_describing_exon_in_aligned_seq(aligned_seq)
		end
		@introns.each do |intron|
			intron.set_variables_describing_intron_in_aligned_seq(aligned_seq)
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

	def get_all_gaps_in_aligned_seq_within_range
		all_gap_pos_with_length = {}

		get_aligned_seq_within_range.each_char.with_index do |chr, ind|
			if chr == "-" then 
				# its a gap

				# special case: the sequence starts with an gap
				if ind == 0 then 
					all_gap_pos_with_length[ ind ] = 1
					next
				end

				last_chr = @aligned_seq[ind - 1]
				if last_chr == "-" then 
					# gap extension
					last_gap_pos = all_gap_pos_with_length.keys.max
					all_gap_pos_with_length[ last_gap_pos ] += 1
				else
					# gap start
					all_gap_pos_with_length[ ind ] = 1
				end
			end
		end
		# list is already reduced to alignment range, as it was iterated over reduced alignment
		return all_gap_pos_with_length
	end

	def get_all_exons_with_length(is_convert_to_nt_length=false)
		# collect start and lenght of each exon
		if is_convert_to_nt_length then 
			all_pos_with_length = @exons.collect do |exon|
				[exon.start_pos_in_dna_seq, exon.length_in_nt]
			end
		else
			all_pos_with_length = @exons.collect do |exon|
				[exon.start_pos_in_aligned_protein, exon.length_in_alignment]
			end 
		end
		return reduce_list_to_range( all_pos_with_length, is_convert_to_nt_length=false )
	end

	def get_all_introns_with_length(is_convert_to_nt_length=false)
		if is_convert_to_nt_length then
			all_pos_with_length = @introns.collect do |intron|
				[intron.pos_last_nt_in_dna_seq_before_intron, intron.n_nucleotides]
			end
		else
			all_pos_with_length = @introns.collect do |intron|
				[intron.pos_last_aa_in_aligned_protein_before_intron, intron.n_nucleotides]
			end
		end
		return reduce_list_to_range( all_pos_with_length, is_convert_to_nt_length=false )
	end

	def get_all_introns_with_phase

		all_pos_with_phase = @introns.collect do |intron|
			intron.get_alignmentpos_and_phase
		end
		return reduce_list_to_range( all_pos_with_phase )
	end

	# # a conserved intron is at same position and phase as an intron in another gene
	# def get_all_conserved_introns
	# 	@introns.select do |intron|
	# 		intron.is_conserved
	# 	end
	# end

	# a sequence of same length as "@aligned_seq" consisting of gaps and intron phases only
	# exon_representation will be used to display exon, default: "-"
	# intron_representation will be used to display intron, default: intron phase
	def plot_intron_phases_onto_aligned_seq(exon_representation, intron_representation)
		exon_representation ||= "-"
		intron_representation ||= nil

		intron_pos_in_alignment = Array.new(@aligned_seq.size, exon_representation)
		@introns.each do |intron|
			pos = intron.pos_last_aa_in_aligned_protein_before_intron
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

	# reduce a list of positions and features to @alignment_range
	# return everything of list which position is in an range
	def reduce_list_to_range(list, is_pos_converted_to_nt=false)
		reduced_list = []
		@alignment_range.each do |start_at, stop_at|
			reduced_list << list.select do |pos_feature|
				start_at <= pos_feature[0] && pos_feature[0] < stop_at
			end
		end
		return reduced_list
	end

end