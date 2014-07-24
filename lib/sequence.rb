module Sequence
	extend self

	# optional parameter: is_return_pos_before_gap: return position of last aa before gap, instead of last gap pos if applicable
	def sequence_pos2alignment_pos(spos, aseq, is_return_pos_before_gap=true)
		pats = []
		aseq.gsub("-", "")[0..spos].split("").each {|chr| pats << ("-*" + chr)}
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
		used_seqs.map do |seq|
			this_length = seq.size
			if this_length < max_length then
				# extend seq with gaps
				seq = seq.ljust(max_length, "-")
			end
		end

		return used_seqs
	end

	# remove common gaps from seqs - array
	# seqs [Array]: array of strings, each string must be of same length
	# opts [Hash]: options, each has a default value
	# 			key gap_symbol [String]: symbol treated as common gap. Must have length 1, default: "-"
	# 			key start_col [Integer]: First column to be searched for common gaps, default: 0
	# 			key delete_all_common_gaps [Boolean]: delete all common gaps or keep one of each consecutive common gaps, default: false
	# 			key ensure_common_gap_between_consecutive_non_gaps [Boolean]: always have one common gap between two non-gaps, default: true
	# 																			only applicable in combination with delete_all_common_gaps == false
	def remove_common_gaps(seqs, opts={})
		gap_symbol = opts[:gap_symbol] || "-"
		start_col = opts[:start_col] || 0
		is_delete_all_common_gaps = false # default: keep one common gap of consecutive ones
		is_always_common_gap_between_two_non_gaps = nil # this option becomes only active when is_delete_all_common_gaps = false
		if opts[:delete_all_common_gaps] then 
			is_delete_all_common_gaps = true
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

		(last_col-1).downto(start_col) do |col|
			content_this_col = reduced_seqs.collect { |seq| seq[col] }
			is_only_gaps_this_col = is_array_of_common_gaps(content_this_col, gap_symbol)

			if is_delete_all_common_gaps then 
				# delete all columns of common gaps, no matter what

				if is_only_gaps_this_col then 
					reduced_seqs.each { |seq| seq.slice!(col) }
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

		return reduced_seqs
	end

	def is_array_of_common_gaps(arr, gap_symbol="-")

		uniq_ele = arr.uniq
		# true if uniq_ele contains only one element, which is gap_symbol
		# false otherwise
		return uniq_ele == [gap_symbol]

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
	end

end