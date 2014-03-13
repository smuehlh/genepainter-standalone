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
	
	# remove common gaps from seqs-array
	# input:
	# seqs [Array] Array of strings, strings should be of same length!
	# opts [Hash]: keep_one_common_gap_of_each_set [Boolean] -> default true
	# 	if true, one common gap of a sequence of common gaps is left and gaps are created whereever necessary
	# 	if false, all common gaps are removed
	# opts [Hash]: insert_common_gap_between_non_gaps [Boolean] -> defaults to keep_one_common_gap_of_each_set
	# 	if true, one common gap is inserted between each consecutive cols being no common-gap cols
	# 	if false, no such gaps are inserted
	# opts [Hash]: gap_symbol [String] -> default: "-" The symbol treated as common gap
	# opts [Hash]: start_col [Integer] -> default: 0 The first column to be searched for common gaps
	def remove_common_gaps(seqs, opts = {})

		# set options or default parameters
		if opts[:keep_one_common_gap_of_each_set].nil? then
			is_keep_one_gap = true
		else
			is_keep_one_gap = opts[:keep_one_common_gap_of_each_set]
		end
		if opts[:insert_common_gap_between_non_gaps].nil? then
			is_add_one_gap = is_keep_one_gap
		else
			is_add_one_gap = opts[:insert_common_gap_between_non_gaps]
		end
		gap_symbol = opts[:gap_symbol] || "-"
		start_col = opts[:start_col] || 0


		last_col = seqs.first.size # all seqs should be of same size ...

		is_only_gaps_this_col, is_only_gaps_after_this_col = nil, nil # "col after"

		# inspect ever pair of consecutive columns
		# find out if current col contains only gaps, and delete it if so
		(last_col-1).downto(start_col) do |col|

			content_this_col = seqs.collect { |seq| seq[col] }
			is_only_gaps_this_col = is_array_of_common_gaps(content_this_col, gap_symbol)

			if is_only_gaps_this_col then 
				# column contains only gaps. check requirements if it should be deleted

				if is_keep_one_gap then 
					#  delete only those columns of common gaps, which are followed by common gap 
					if is_only_gaps_after_this_col then 
						seqs.each { |seq| seq.slice!(col) }
					end
				else
					# delete every column of common gaps, regardless the surronding columns
					seqs.each { |seq| seq.slice!(col) }
				end
			else
				# column contains also non-gaps. check if a column of common gaps should be added
				if is_add_one_gap && ! is_only_gaps_after_this_col then 
					seqs.each { |seq| seq.insert(col+1, gap_symbol) }	
				end
			end

			# this column will be column after the current column in next iteration
			is_only_gaps_after_this_col = is_only_gaps_this_col
		end

		return seqs
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