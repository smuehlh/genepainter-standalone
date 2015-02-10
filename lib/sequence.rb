module Sequence
	extend self

 	# optional parameter: is_return_pos_before_gap: return position of last aa before gap, instead of last gap pos if applicable
    def sequence_pos2alignment_pos(spos, aseq, is_return_pos_before_gap=true)
        pats = []

        aseq_wo_gaps = aseq.delete("-")
		aseq_wo_gaps[0..spos].split("").each do |chr| 
			pat = "-*"
			if chr == "*" then 
				# a stop codon (denoted by '*'), escape it!
				pat += "\\*"
			else
				pat += chr
			end
			pats << pat
		end
       	if ! is_return_pos_before_gap then 
                # add pattern for trailing gap
                pats.push("-*")
        end
        pat = Regexp.new(pats.join)
        pat.match(aseq)[0].length - 1
    end

	def alignment_pos2sequence_pos(apos, aseq)
		mapped_pos = aseq[0..apos].gsub("-", "").length - 1
		# gaps at the begin of the sequence are mapped onto position - 1
		# shift them back to position 0
		mapped_pos = 0 if mapped_pos < 0
		return mapped_pos
	end

	def is_protein_sequence_valid(seq)
		if seq.match(/[^A-Z*-]/i) || seq.empty? then
			# sequence contains special chars other than "-" (gap symbol), "*" (stop codon)
			return false
		else
			return true
		end
	end	
	def replace_invalid_chars_in_protein_sequence(seq)
		return seq.gsub(/[^A-Z*-]/i, "X")
	end

	# read in alignment, splits file into fasta headers and sequences
	# removes any species names from fasta headers and returns them separately
	def read_in_alignment(path)
		names, seqs = [], []
		IO.foreach(path) do |line|
			line.chomp!
			if line.start_with? ">" then
				# fasta header
				names << line[1..-1]
			else
				# fasta sequence
				n_seqs = names.size
				# a new sequence or another line for the last sequence?
				if seqs.size < n_seqs then
					# new sequence
					seqs << line
				else
					# add to last sequence
					seqs[n_seqs - 1] += line # -1: convert number of elements n_seqs to an index
				end
			end
		end
		names.delete("")
		seqs.delete("")
		if names.size != seqs.size then
			Helper.abort "Error while parsing multiple sequence alignment. Number of fasta header and sequences does not match." 
		end
		return names, seqs

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "read_in_alignment", exp
    	Helper.abort "Cannot read alignment."
	end

	def convert_strings_to_fasta(fasta_header, fasta_seq)
		fasta_arr = []

		# if necessary, add ">" at beginning of header to make it valid
		fasta_header = ">" << fasta_header if ! fasta_header.start_with?(">")
		fasta_arr << fasta_header

		fasta_arr << fasta_seq.scan(/.{1,80}/)

		return fasta_arr.join("\n")
	end

	# ensure that all sequences, which will be finally used (!) - and only those - have same lenght
	def ensure_seqs_have_same_length(aligned_seqs, seqs_names, used_names)
		# find max length of all sequences which will be finally used
		max_length = 0
		used_seqs = []

		# used_seqs is of same order as used_names
		used_names.each do |name|
			seqs_ind = seqs_names.index(name)
			if seqs_ind then
				seq = aligned_seqs[seqs_ind] # aligned_seqs is of same order as seqs_names
				used_seqs << seq
				this_length = seq.size
				if max_length < this_length then
					max_length = this_length
				end
			end
		end

		# make sure every seq, which is finally used, is of same length
		used_seqs.map! do |seq|
			this_length = seq.size
			if this_length < max_length then
				# extend seq with gaps
				seq.ljust(max_length, "-")
			else
				seq
			end
		end

		return used_seqs
    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "ensure_seqs_have_same_length", exp
    	Helper.abort "Cannot ensure all aligned sequences have same length."		
	end

	# seqs [Array]: array of strings, each must be of same length
	# range [Hash]: range parameters
		# key reverser_position [Array]: List of positions in sequence
	# range[:reverse_position] might contain keyword Infinity, which stands for "last pos in alignment"
	# replace this keyword by the actuall last position
	def replace_range_keyword_end_of_alignment_by_position(seqs, range)
		if range[:reverse_position] && range[:reverse_position][0] then 
			range[:reverse_position][0] = seqs[0].size - 1 # -1 as seqs[seqs.size] is nil
		end
	end

	def remove_common_gaps_in_alignment_update_predefined_ranges(seqs, range)
		reduced_seqs, deleted_pos = remove_common_gaps(seqs, 
			{
				intron_pattern: false, # its an protein alignment, no exon-intron pattern
				delete_all_common_gaps: true, # remove all common gaps
				return_deleted_cols: true # return indices of columns that were deleted
			}
			)

		range_pos_in_reduced_seq = []
		range[:reverse_position].each_slice(2) do |range_end, range_start|
			n_deleted_cols_until_range_start = deleted_pos.select{|pos| pos < range_start}.size
			n_deleted_cols_until_range_end = deleted_pos.select{|pos| pos <= range_end}.size

			range_pos_in_reduced_seq.push range_end - n_deleted_cols_until_range_end
			range_pos_in_reduced_seq.push range_start - n_deleted_cols_until_range_start
		end

		range[:reverse_position] = range_pos_in_reduced_seq

		return reduced_seqs
	end

	# remove common gaps from seqs - array
	# seqs [Array]: array of strings, each string must be of same length
	# opts [Hash]: options, each has a default value
	# 			key gap_symbol [String]: symbol treated as common gap. Must have length 1, default: "-"
	# 			key start_col [Integer]: First column to be searched for common gaps, default: 0
	# 			key intron_pattern [Boolean]: treat input as intron_pattern, 
	# 											i.e. different non-gap symbols in same row are distributed onto several rows, default: true
	# 			key delete_all_common_gaps [Boolean]: delete all common gaps or keep one of each consecutive common gaps, default: false
	# 			key ensure_common_gap_between_consecutive_non_gaps [Boolean]: always have one common gap between two non-gaps, default: true
	# 																			only applicable in combination with delete_all_common_gaps == false
	# 			key return_deleted_cols [Boolean]: return array with all deleted columns, default: false 
	# 												only applicable in combination with delete_all_common_gaps == true
	def remove_common_gaps(seqs, opts={})
		gap_symbol = opts[:gap_symbol] || "-"
		start_col = opts[:start_col] || 0
		is_delete_all_common_gaps = false # default: keep one common gap of consecutive ones
		is_always_common_gap_between_two_non_gaps = nil # this option becomes only active when is_delete_all_common_gaps = false
		is_return_deleted_pos = false # this option becomes only active when is_delete_all_common_gaps = true
		if opts[:intron_pattern].nil? then 
			# split rows with different intron phases into different rows
			is_intron_pattern = true
		else
			# _not_ split rows with intron phases into different rows
			is_intron_pattern = opts[:intron_pattern]
		end
		if opts[:delete_all_common_gaps] then 
			is_delete_all_common_gaps = true
			if opts.has_key?(:return_deleted_cols) &&
				opts[:return_deleted_cols] == true then 
				is_return_deleted_pos = true
			end
		else
			is_delete_all_common_gaps = false
			is_always_common_gap_between_two_non_gaps = true # default: always have a common gap between two non-gaps (insert one if need is)
			if opts.has_key?(:ensure_common_gap_between_consecutive_non_gaps) &&
				opts[:ensure_common_gap_between_consecutive_non_gaps] == false then 
				is_always_common_gap_between_two_non_gaps = false
			end
		end

		reduced_seqs = seqs.map {|ele| ele.dup} # iterate of this set of sequences to keep the input seqs unchanged
		last_col = reduced_seqs.first.size
		last_col_to_keep = last_col-1 # -1: index starts with 0
		deleted_pos = []

		(last_col-1).downto(start_col) do |col|
			content_this_col = reduced_seqs.collect { |seq| seq[col] }
			is_only_gaps_this_col = is_array_of_common_gaps(content_this_col, gap_symbol)
			intron_phases_this_col = intron_phases_in_array(content_this_col)
			is_multiple_intron_phases_this_col = ! is_only_gaps_this_col && intron_phases_this_col.size > 1

			if is_delete_all_common_gaps then 
				# delete all columns of common gaps, no matter what

				if is_only_gaps_this_col then 
					reduced_seqs.each { |seq| seq.slice!(col) }
					deleted_pos.push col
				end
			else
				# delete all consecutive common gaps but one
			
				if ! is_only_gaps_this_col then
					# column contains also non-gaps

					# reduce consecutive common-gap cols to one colum
					n_cols_to_del = last_col_to_keep - col - 1 # -1: keep always one col
					reduced_seqs.each{ |seq| seq.slice!(col+1,n_cols_to_del) }

					if is_always_common_gap_between_two_non_gaps && 
						( col == last_col_to_keep ) then 

						# this is a non-gap col, after this col in string is another non-gap col: separate them by common gap
						reduced_seqs.each { |seq| seq.insert(col+1, gap_symbol) }
					end
					if ! is_always_common_gap_between_two_non_gaps then 
	
						# this is a non-gap col, 2 cols after this col (in string) is another non-gap col
						# if the non-gap is not in same sequence, delete the common-gap-col in between
						is_non_gap_in_same_seq = false
						reduced_seqs.each do |seq|
							if seq[col] != gap_symbol && seq[col+2] && seq[col+2] != gap_symbol then 
								is_non_gap_in_same_seq = true
								break
							end
						end
						if ! is_non_gap_in_same_seq then 
							# its safe to delete common gap in between
							reduced_seqs.each{ |seq| seq.slice!(col+1) }
						end
					end

					# set the col to keep to column preceeding (in the string) this one!
					last_col_to_keep = col - 1 
				end
			end	

			if is_intron_pattern && is_multiple_intron_phases_this_col then 
				# column contains different intron phases
				# separate different intron phases onto different columns
				# insert new columns behind this one

				# intron_phases_this_col are sorted by intron phase (= numerically sorted)
				n_intron_phases = intron_phases_this_col.size
				if is_delete_all_common_gaps then 
					# split introns onto different cols, no need to add common gap cols in between
					separator = gap_symbol * n_intron_phases
				else
					# split introns onto different cols, add common gap cols in between
					separator = gap_symbol * (n_intron_phases * 2 - 1) # -1: col for first intron already exists; *2 for common gaps col
				end
				reduced_seqs.each do |seq|
					this_char = seq[col]
					this_separator = separator.dup # preserve original separator
					if this_char != gap_symbol then 
						# insert actual intron phase in separator
						ind_intron_phase = intron_phases_this_col.index(this_char)
						if ! is_delete_all_common_gaps then 
							ind_intron_phase *= 2 # account for common gap cols between each two non-gap cols
						end
						this_separator[ind_intron_phase] = this_char
					end

					# replace col by separator
					seq[col] = this_separator
				end
			end

		end
		if ! is_delete_all_common_gaps then 
			# reduce consecutive common gaps to one col
			n_cols_to_del = last_col_to_keep - start_col # not -1 because col itself it the column to always keep
			reduced_seqs.each { |seq| seq.slice!(start_col,n_cols_to_del) }

			# ensure that reduced_seqs always start with common gap
			first_col = reduced_seqs.collect { |seq| seq[start_col] }
			if ! is_array_of_common_gaps(first_col, gap_symbol) then 
				reduced_seqs.each { |seq| seq.insert(start_col, gap_symbol) }
			end

			# ensure that reduced_seqs always end with common gap
			last_col = reduced_seqs.collect { |seq| seq[-1] }
			if ! is_array_of_common_gaps(last_col, gap_symbol) then 
				reduced_seqs.each { |seq| seq.insert(-1, gap_symbol) }
			end
		end

		if is_return_deleted_pos then 
			return reduced_seqs, deleted_pos
		else
			return reduced_seqs
		end
    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "remove_common_gaps", exp
    	Helper.abort "Cannot remove common gaps in alignment."

	end

	def is_array_of_common_gaps(arr, gap_symbol="-")

		uniq_ele = arr.uniq
		# true if uniq_ele contains only one element, which is gap_symbol
		# false otherwise
		return uniq_ele == [gap_symbol]
	end

	def intron_phases_in_array(arr)
		# find all intron phases in array
		return arr & ["0","1","2","?"]
	end

	# map position of intron onto reduced pattern
	# do not use intron_col * 2 + 1, since reduced pattern might not always contain a common-gap col between two intron cols
	def find_index_of_intron_pos(pattern, intron_col)
		n_intron_cols_seen = 0
		n_cols = pattern[0].size
		0.upto(n_cols-1) do |col|
			content_this_col = pattern.collect { |seq| seq[col] }
			is_only_gaps_this_col = is_array_of_common_gaps(content_this_col)
		
			if ! is_only_gaps_this_col then 
				return col if n_intron_cols_seen == intron_col

				n_intron_cols_seen += 1
			end
		end
	end

	# check if there is at least one sequence without gap at this position
	# check if that sequence has an intron at gap start/end?
	# re-position intron for each gene in genes-lists
	def correct_introns_flanking_gaps(introns_before_gap_pos_gene, gene_objects)

		introns_before_gap_pos_gene.each do |pos_before_gap, occurence|

			gene_names_this_intron = occurence.collect { |arr| arr[0] }
			gap_end_positions_this_intron = occurence.collect { |arr| arr[1] }

			gene_objects.each do |ref_gene|

				if gene_names_this_intron.include?(ref_gene.name) then 
					# reference gene has a gap itself after this intron
					next
				end

				# reference gene has no gaps starting at position of this intron

				# introns of this gene might be at same position as intron under question:
				# 1) at pos_before_gap: don't do anything
				# 2) at pos after gap: change position of introns!

				intron_pos_ref_gene = ref_gene.get_all_intronpositions

				gap_end_positions_this_intron.each_with_index do |pos_end_of_gap, ind|

					if intron_pos_ref_gene.include?(pos_end_of_gap) then 

						# shift intron of this gene to end of gap
						this_name = gene_names_this_intron[ind]
						this_gene = gene_objects.select { |g| g.name ==  this_name }.first # select returns array, use first to get object itself

						this_gene.introns.each do |intron|
							if intron.pos_last_aa_in_aligned_protein_before_intron == pos_before_gap then 
								# this is the intron under question
								intron.set_variables_describing_intron_in_aligned_seq( this_gene.aligned_seq, false ) #false: do not position intron at beginning of gap
							end
						end
					end
				end
			end
		end

    rescue NoMethodError, TypeError, NameError, ArgumentError, Errno::ENOENT => exp
    	Helper.log_error "correct_introns_flanking_gaps", exp
    	Helper.abort "Cannot correct intron positions flanking alignment gaps."		
	end

end