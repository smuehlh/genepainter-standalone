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
#		@alignment_range = [] # array of arrays; restrict mapping of gene stucture onto alignment to these parts of the alignment
	end

	# set instance variable @aligned_seq and also exon/intron position in alignment
	def add_aligned_seq(aligned_seq)
		@aligned_seq = aligned_seq
		@exons.each do |exon|
			exon.set_variables_describing_exon_in_aligned_seq(aligned_seq)
		end
		@introns.each do |intron|
			intron.set_variables_describing_intron_in_aligned_seq(aligned_seq)
		end
	end

	# range might be positive or negative range: keep only range or keep everything but range
	def reduce_gene_to_range(range)
		range_start = range[:position][0]
		range_end = range[:position][1]

		# sanity check: range end must not be outside aligned sequence
		if range_end >= @aligned_seq.size then 
			Helper.abort "end position of range exceeds alignment."
		end

		introns_within_range = []
		exons_within_range = []
		aligned_seq_within_range = ""

		if range[:is_delete_range] then 

			# maybe variable name is abit missleading, but aligned_seq_within_range is sequence _outside range to del_
			aligned_seq_within_range = @aligned_seq[0...range_start] + @aligned_seq[range_end+1..-1]
			aligned_seq_deleted = @aligned_seq[range_start..range_end] # inclusive range_end

			n_del_nt = aligned_seq_deleted.delete("-").size * 3

			is_continue_exon_after_range = false # kept exon ending in range and kept exon starting in range must be merged!

			@exons.each do |exon|

				# exon is within parts to keep when:
				# start pos is before range (to delete); end pos: before/in/after
				# OR end pos is after range (to delete); start pos: in/after

				is_start_pos_before_range = exon.start_pos_in_aligned_protein <= range_start
				is_start_pos_within_range = ( range_start <= exon.start_pos_in_aligned_protein && exon.start_pos_in_aligned_protein <= range_end)

				is_end_pos_before_range = exon.end_pos_in_aligned_protein <= range_start
				is_end_pos_within_range = ( range_start <= exon.end_pos_in_aligned_protein && exon.end_pos_in_aligned_protein <= range_end )

				starting_somewhere_special = nil
				ending_somewhere_special = nil

				is_collect_exon = false

				if is_start_pos_before_range then 

					# keep start pos
					starting_somewhere_special = exon.start_pos_in_dna_seq

					if is_end_pos_before_range then 

						# change nothing
						ending_somewhere_special = exon.end_pos_in_dna_seq
						is_collect_exon = true

					elsif is_end_pos_within_range 

						ending_somewhere_special = @aligned_seq[0..range_start].delete("-").size * 3
						is_collect_exon = true
						is_continue_exon_after_range = true

					else # is_end_pos_after_range
					
						is_collect_exon = true

					end

				elsif is_start_pos_within_range

					# change start and end pos

					if is_end_pos_within_range then 

						# nothing to do; don't keep exon.

					else # is_end_pos_after_range

						if is_continue_exon_after_range then 
							last_kept_exon = exons_within_range.pop

							starting_somewhere_special = last_kept_exon.start_pos_in_dna_seq
							is_continue_exon_after_range = false
						else
							starting_somewhere_special = @aligned_seq[0..range_start].delete("-").size * 3 
						end

						is_collect_exon = true

					end

				else # is_start_pos_after_range

					# is_end_after_range

					is_collect_exon = true

				end
					
							
				if is_collect_exon then 
					exon_conform_to_range = exon.create_copy_with_shifted_positions(n_del_nt, 
						starting_somewhere_special, ending_somewhere_special, 
						aligned_seq_within_range)
					exons_within_range.push( exon_conform_to_range )
				end						

			end # exons.each

			# intron is within parts to keep when:
			# pos last aa before intron is before range or after range
			# collect introns

			@introns.each do |intron|

				if intron.pos_last_aa_in_aligned_protein_before_intron <= range_start then 
					# intron is before range

					# do not change intron position
					introns_within_range.push( intron )

				elsif intron.pos_last_aa_in_aligned_protein_before_intron >= range_end

					intron_conform_to_range = intron.create_copy_with_shifted_positions( n_del_nt, aligned_seq_within_range) 
					introns_within_range.push( intron_conform_to_range )

				else
					# intron is within range, nothing to do

				end

			end # introns.each

		else

			aligned_seq_within_range =  @aligned_seq[range_start..range_end] # inclusive range_end
			aligned_seq_before_range = @aligned_seq[0...range_start]
			n_del_nt = aligned_seq_before_range.delete("-").size * 3 # exclusive range_start


			# collect exons and introns within range
			@exons.each do |exon|

				# exon is within range when:
				# start pos of exon is within range
				# OR end pos of exon is within range
				# OR both is within range

				is_start_pos_within_range = ( range_start <= exon.start_pos_in_aligned_protein && exon.start_pos_in_aligned_protein <= range_end)
				is_end_pos_within_range = ( range_start <= exon.end_pos_in_aligned_protein && exon.end_pos_in_aligned_protein <= range_end )

				starting_somewhere_special = nil # needed if exon starts outside range
				ending_somewhere_special = nil # needed if exon ends outside range

				is_collect_exon = false

				if is_start_pos_within_range && is_end_pos_within_range then 
					is_collect_exon = true

				elsif is_start_pos_within_range 
					# end is outside range; set it to range border
					ending_somewhere_special = aligned_seq_within_range.delete("-").size * 3
					is_collect_exon = true

				elsif is_end_pos_within_range
					# start is outside range; set it to range border
					starting_somewhere_special = 0
					is_collect_exon = true

				elsif exon.start_pos_in_aligned_protein <= range_start && range_end <= exon.end_pos_in_aligned_protein
					# exon spans the range
					# start outside range; set start to border
					starting_somewhere_special = 0
					# ends outside range; set end to border
					ending_somewhere_special = aligned_seq_within_range.delete("-").size * 3
					is_collect_exon = true
					
				else
					# exon is not in range
				end		

				if is_collect_exon then 
					exon_conform_to_range = exon.create_copy_with_shifted_positions(n_del_nt, 
						starting_somewhere_special, ending_somewhere_special, 
						aligned_seq_within_range)
					exons_within_range.push( exon_conform_to_range )
				end

			end # exons.each


			# intron is within range when:
			# pos last aa before intron is within range
			# collect introns

			@introns.each do |intron|

				if range_start <= intron.pos_last_aa_in_aligned_protein_before_intron &&
					intron.pos_last_aa_in_aligned_protein_before_intron <= range_end then 

					intron_conform_to_range = intron.create_copy_with_shifted_positions(n_del_nt, aligned_seq_within_range)
					introns_within_range.push( intron_conform_to_range )
				end

			end # introns.each

		end

		@aligned_seq = aligned_seq_within_range
		@exons = exons_within_range
		@introns = introns_within_range


		if @aligned_seq.empty? ||
			@exons.empty? then  
			# allow 1-exon genes (no introns), but expects at least 1-exon genes
			Helper.abort "Cannot reduce gene #{@name} to range."
		end
	end


	# returns common introns of both genes, thus, returned introns are not uniq
	def common_introns_of_this_and_other_gene(other_gene)
		(@introns & other_gene.introns).concat(other_gene.introns & @introns)
	end

	def get_all_gaps_in_aligned_seq
		all_gap_pos_with_length = {}

		@aligned_seq.each_char.with_index do |chr, ind|
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
		return all_pos_with_length
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
		return all_pos_with_length
	end

	def get_all_introns_with_phase

		all_pos_with_phase = @introns.collect do |intron|
			intron.get_alignmentpos_and_phase
		end
		return all_pos_with_phase
	end

	def get_all_intronpositions
		@introns.collect do |intron|
			intron.pos_last_aa_in_aligned_protein_before_intron
		end
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

		reduced_pattern = intron_pos_in_alignment.join("")
		return reduced_pattern
	end

end