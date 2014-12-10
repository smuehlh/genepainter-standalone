# a gene consists of exons and introns

class Gene
	attr_accessor :aligned_seq, :exons, :introns, :taxonomic_lineage, :ind_first_uniq_ancestor
	attr_reader :name

	def initialize(name)
		@name = name
		@aligned_seq = nil
		@exons = [] # exon objects in their correct order
		@introns = [] # intron objects in their correct order

		@taxonomic_lineage = []
		@ind_first_uniq_ancestor = nil
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

	# lineage: array, starting with root
	# ind: the index of the first uniq ancestor of the species in that array -> later replaced by the array element
	def add_taxonomy(lineage_and_ind_first_uniq_ancestor)
		lineage = lineage_and_ind_first_uniq_ancestor[:lineage]
		ind = lineage_and_ind_first_uniq_ancestor[:ind_first_uniq]
		@taxonomic_lineage = lineage
		@ind_first_uniq_ancestor = ind 
	end

	def reduce_gene_to_range(range)
		# range_end might be Infinity, which means the last aligment position should be used!
		is_delete_range = range[:is_delete_range]

		range[:reverse_position].each_slice(2) do |range_end, range_start|
			reduce_gene_to_single_range( range_start, range_end, is_delete_range )
		end
	end

	# range might be positive or negative range: keep only range or keep everything but range
	def reduce_gene_to_single_range(range_start, range_end, is_delete_range)

		# sanity check: range end must not be outside aligned sequence
		if range_end >= @aligned_seq.size then 
			Helper.abort "end position of range #{range_start}-#{range_end} exceeds alignment."
		end

		introns_within_range = []
		exons_within_range = []
		aligned_seq_within_range = ""

		if is_delete_range then 
			# delete range

			# maybe variable name is a bit missleading, but aligned_seq_within_range is sequence _outside range to del_
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
			# keep range

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
puts "#{range_start}-#{range_end}"
puts @aligned_seq
		if @aligned_seq.empty? ||
			@exons.empty? then  
			# every gene should have still min. 1 exon
			Helper.warn "Cannot reduce gene #{@name} to range."
			throw :error
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

	def get_all_exons_with_length
		# collect start and lenght of each exon
		@exons.collect do |exon|
			[exon.start_pos_in_aligned_protein, exon.length_in_alignment]
		end 
	end

	def get_all_introns_with_length(is_convert_to_nt_length=false)
		@introns.collect do |intron|
			[intron.pos_last_aa_in_aligned_protein_before_intron, intron.n_nucleotides]
		end
	end

	def get_all_introns_with_phase

		@introns.collect do |intron|
			intron.get_alignmentpos_and_phase
		end
	end

	def get_all_intronpositions
		@introns.collect do |intron|
			intron.pos_last_aa_in_aligned_protein_before_intron
		end
	end

	def get_all_intronpositions_merged_with_phase
		get_all_introns_with_phase.collect do |pos, phase|
			Intron.merge_position_and_phase(pos, phase)
		end
	end

	# # a conserved intron is at same position and phase as an intron in another gene
	# def get_all_conserved_introns
	# 	@introns.select do |intron|
	# 		intron.is_conserved
	# 	end
	# end

	def get_all_gap_boundaries_preceeded_by_intron
		@introns.collect do |intron|

			char_after_intron_in_alignment = @aligned_seq[ intron.pos_last_aa_in_aligned_protein_before_intron + 1]
			if char_after_intron_in_alignment == "-" then
				# intron is located before an gap

				# but due to some very strange gene prediction, the intron might split the very last codon (yes, seen once!)
				begin
					gap = @aligned_seq[ intron.pos_last_aa_in_aligned_protein_before_intron + 1 .. -1].match(/(-+)[^-]/)[1]
				rescue
					gap = ""
				end
				pos_gap_end = intron.pos_last_aa_in_aligned_protein_before_intron + gap.size # last position of gap

				[intron.pos_last_aa_in_aligned_protein_before_intron, pos_gap_end]
			else
				nil
			end
		end.compact
	end

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

		intron_pos_in_alignment.join("")
	end

	def last_common_ancestor_with_other_gene(other_lineage)
		# as lineages start with root, the last common element is the last common ancestor
		return common_lineage_with_other_gene(other_lineage).last
	end

	def common_lineage_with_other_gene(other_lineage)
		return @taxonomic_lineage.intersection( other_lineage )
	end
	
	def get_lineage_root_to_first_uniq_ancestor_of_species
		if has_taxonomic_information then 
			return @taxonomic_lineage[0..@ind_first_uniq_ancestor] 
		else
			return []
		end
	end

	def first_uniq_ancestor_of_species
		if has_taxonomic_information then 
			return @taxonomic_lineage[@ind_first_uniq_ancestor]
		else
			return ""
		end
	end
	def has_taxonomic_information
		@taxonomic_lineage.any? && @ind_first_uniq_ancestor
	end

end