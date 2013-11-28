module Sequence
	extend self

	def sequence_pos2alignment_pos(spos, aseq)
		pats = []
		aseq.gsub("-", "")[0..spos].split("").each {|chr| pats << ("-*" + chr)}
		pat = Regexp.new(pats.join)
		pat.match(aseq)[0].length - 1
	end

	def alignment_pos2sequence_pos(apos, aseq)
		aseq[0..apos].gsub("-", "").length - 1
	end

	def convert_strings_to_fasta(fasta_header, fasta_seq)
		fasta_arr = []

		# if necessary, add ">" at beginning of header to make it valid
		fasta_header = ">" << fasta_header if ! fasta_header.start_with?(">")
		fasta_arr << fasta_header

		fasta_arr << fasta_seq.scan(/.{1,80}/)

		return fasta_arr.join("\n")
	end

	# # input: array containing aligned sequences
	# def remove_common_gaps_from_aligned_seqs(aligned_seqs)

	# 	max_length = aligned_seqs.collect { |s| s.size }.max 
	# 	n_aligned_seqs = aligned_seqs.size

	# 	# look for common gaps
	# 	gaps = []
	# 	max_length.times do |pos|
	# 		n_seqs_with_gap = aligned_seqs.select { |seq| seq[pos] && seq[pos] == "-" }.length
	# 		if n_seqs_with_gap == n_aligned_seqs then
	# 			gaps << pos
	# 		end
	# 	end
	# 	gaps = gaps.reverse # otherwise delete_at will not work as expected

	# 	# delete common gaps
	# 	aligned_seqs_without_common_gaps = aligned_seqs.collect do |seq|
	# 		seq_splitted_into_chars = seq.split("")
	# 		gaps.each { |pos| seq_splitted_into_chars.delete_at(pos) }
	# 		seq_splitted_into_chars.join("")
	# 	end

	# 	return aligned_seqs_without_common_gaps

	# end


	# remove not all common gaps but only those which are not needed for gene structure 
	# they are not needed, if next column also contains only gaps
	def remove_common_gaps_from_gene_structurs(exon_intron_patterns)

		max_length = exon_intron_patterns.collect { |s| s.size }.max 

	# TODO will be faster by thinking of exon_intron_patterns as an array of arrays
	# & representing it as hash or matrix !!!
	# then, it might be possible to do both in one step

		# look for common gaps which are not neccessary for the exon_intron_pattern
		gaps = []
		(max_length).times do |pos|

			# this will fail if the size of the array is exceeded ???
			pattern_at_column = exon_intron_patterns.collect { |seq| seq[pos] }
			pattern_at_next_column = exon_intron_patterns.collect { |seq| seq[pos+1] }

			if is_array_of_common_gaps(pattern_at_column) then
				# column contains only gaps
				if is_array_of_common_gaps(pattern_at_next_column) then
					# next column also contains only gaps
					# its save to remove column

					gaps.unshift(pos)
				end 
			end
		end

		# delete common gaps
		reduced_exon_intron_patterns = exon_intron_patterns.collect do |seq|
			seq_splitted_into_chars = seq.split("")
			gaps.each { |pos| seq_splitted_into_chars.delete_at(pos) }
			seq_splitted_into_chars.join("")
		end		

		return reduced_exon_intron_patterns

	end

	def is_array_of_common_gaps(arr)

		uniq_ele = arr.uniq
		if uniq_ele.size == 1 && uniq_ele.include?("-") then
			return true
		end
		return false
	end

end