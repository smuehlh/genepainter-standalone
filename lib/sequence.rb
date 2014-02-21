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


	# remove common gaps from seqs-array, after position offset
	# input "seqs" can be alignment or gene structures
	# input "seqs" _must be_ of same length
	# in case of alignment, all common gaps are removed
	# in case of exon_intron_pattern, only those gaps are removed, where the next column also contains only common gaps
	# opts [Hash] Boolean is_alignment -> default: false
	# opts [Hash] Integer start_col -> default: 0
	# opts [Hash] String gap_symbol -> default: "-"
	def remove_common_gaps(seqs, opts = {})
		# set options or default params
		is_alignment = opts[:is_alignment] || false
		start_col = opts[:start_col] || 0
		gap_symbol = opts[:gap_symbol] || "-"

		last_col = seqs.collect { |s| s.size }.max

		is_only_gaps_next_col = false
		is_only_gaps_this_col = false

		# inspect ever pair of consecutive columns
		# find out if current col contains only gaps, and delete it if so
		(last_col-1).downto(start_col) do |col|
			content_of_this_col = seqs.collect { |seq| seq[col] }

			is_only_gaps_this_col = is_array_of_common_gaps(content_of_this_col, gap_symbol)

			if is_only_gaps_this_col then
				if is_alignment then
					# delete every column of common gaps, regardless the surronding columns
					seqs.each { |seq| seq.slice!(col) }
				else
					# delete only those columns of common gaps, which are followed by common gap 
					if is_only_gaps_next_col then
						seqs.each { |seq| seq.slice!(col) }
					end
				end
			end

			# this column will be last column in next iteration
			is_only_gaps_next_col = is_only_gaps_this_col
		end

		return seqs
	end

	def is_array_of_common_gaps(arr, gap_symbol="-")

		uniq_ele = arr.uniq
		# true if uniq_ele contains only one element, which is gap_symbol
		# false otherwise
		return uniq_ele == [gap_symbol]

	end

end