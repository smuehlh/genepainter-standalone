class NeedlemanWunsch
	attr_reader :aligned_seq_a, :aligned_seq_b, :score

	def initialize(seq_a, seq_b, penalize_endgaps=false)
		@seq_a = seq_a
		@seq_b = seq_b
@seq_a = "GARFIELDTHECAT"
@seq_b = "GARFIELDLASTCAT"
		@penalize_endgaps = penalize_endgaps || false # penalize end gaps?
		@lin_gap_pen = true # use linear or affin linear gap costs?
		@blosum62 = read_in_blosum62

		@aligned_seq_a, @aligned_seq_b, @score = align
	end

	def align
		scores = build_scoring_matrix
		aligned_a, aligned_b, score = trace_back(scores)
	end


# 	def build_scoring_matrix
# 		# TODO 
# 		# this is actually implementation of gotoh algorithm, so maybe re-name class
# 		m = @seq_a.size - 1
# 		n = @seq_b.size - 1
		
# 		# init
# 		smat = Array.new(m+1) { Array.new(n+1) { 0 } } # m rows, n cols # scores,ends in alignment of seqa, seqb
# 		vmat = Array.new(m+1) { Array.new(n+1) { 0 } } # scores, ends in gap in seq_a
# 		hmat = Array.new(m+1) { Array.new(n+1) { 0 } } # scores, ends in gap in seq_b

# 		(0..n).to_a.each do |j|
# 			vmat[0][j] = -1.0/0.0 # -inf
# 			next if j == 0
# 			hmat[0][j] = gap_func(j)
# 			smat[0][j] = gap_func(j)
# 		end

# 		(0..m).to_a.each do |i|
# 			hmat[i][0] = -1.0/0.0 # -inf
# 			next if i == 0
# 			vmat[i][0] = gap_func(i)
# 			smat[i][0] = gap_func(i)
# 		end

# 		# recursion
# 		(1..m).to_a.each do |i|
# 			(1..n).to_a.each do |j|
# 				hmat[i][j] = [ smat[i][j-1] + gap_start + gap_ext, # begin new insertion
# 					hmat[i][j-1] + gap_ext # continue insertion
# 					].max 
# 				vmat[i][j] = [ smat[i-1][j] + gap_start + gap_ext, # begin new deletion
# 					vmat[i-1][j] + gap_ext # continue deletion
# 					].max
# 				smat[i][j] = [smat[i-1][j-1] + aa_score(@seq_a[i], @seq_b[j]), #match
# 					hmat[i][j], # insertion
# 					vmat[i][j] # deletion
# 					].max
# 			end
# 		end
# 		return smat, vmat, hmat
# 	end
# 	def trace_back(smat, vmat, hmat)
# # TODO does not work
# 		m = @seq_a.size - 1
# 		n = @seq_b.size - 1
# 		aligned_seq_a = []
# 		aligned_seq_b = []

# 		max_score = smat.last.last # lower right corner

# 		last_row = smat.last
# 		last_col = smat.collect {|row| row.last}
# 		i = last_col.each_with_index.max[1] # index of maximal value
# 		j = last_row.each_with_index.max[1]

# 		max_lenght = [m,n].max
# 		aligned_seq_a = [ @seq_a[i+1..m].ljust(max_lenght-m, "-") ]
# 		aligned_seq_b = [ @seq_b[j+1..n].ljust(max_lenght-n, "-") ]
# 		while i > 0 && j > 0 do
# 			case [ smat[i][j], vmat[i][j], hmat[i][j] ].max
# 			when smat[i][j]
# 				# match in i and j
# 				aligned_seq_a.unshift(@seq_a[i])
# 				aligned_seq_b.unshift(@seq_b[j])
# 				i = i-1
# 				j = j-1

# 			when vmat[i][j]
# 				# gap in seq_a
# 				aligned_seq_a.unshift(@seq_a[i])
# 				aligned_seq_b.unshift("-")
# 				i = i-1

# 			when hmat[i][j]
# 				# gap in seq_b
# 				aligned_seq_a.unshift("-")
# 				aligned_seq_b.unshift(@seq_b[j])
# 				j = j-1

# 			end
# 		end

# puts "in trace back"
# 		debugger
# 		puts "bla"
# 	end
# 	def gap_start
# 		return -7
# 	end
# 	def gap_ext
# 		return -1
# 	end
# 	def gap_func(len)
# 		return gap_start + len * gap_ext
# 	end

	def build_scoring_matrix
		n_rows = @seq_b.size
		n_cols = @seq_a.size
		scores = Array.new(n_cols) { Array.new(n_rows) { 0 } } # number rows: seq_b, number cols: seq_a

		# init first row (=all gaps in seq_a)
		(0..n_cols-1).to_a.each do |i|
			scores[i][0] = gap_cost(i + 1 ) # gap at index 0 has length 1, and so on ... 
		end
		# init first col (= all gaps in seq_b)
		(0..n_rows-1).to_a.each do |j|
			scores[0][j] = gap_cost(j + 1 ) # gap at index 0 has lenght 1, and so on ...
		end

		this_gap_len_a, this_gap_len_b = 1, 1 # default: insert 1 gap

# TODO this will leave out first char in seq???
# I und j vertauscht, aber methoden in sich konsistent
		(1..n_cols-1).to_a.each do |i| # iterate over seq a
			(1..n_rows-1).to_a.each do |j| # iterate over seq b
				match = scores[i-1][j-1] + aa_score(@seq_a[i], @seq_b[j])
				del = scores[i-1][j] + gap_cost(this_gap_len_a)
				ins = scores[i][j-1] + gap_cost(this_gap_len_b)
				scores[i][j] = [match, del, ins].max

				if [match, del, ins].max == match then 
					this_gap_len_a, this_gap_len_b = 1, 1
				elsif [match, del, ins].max == del 
					this_gap_len_a += 1
				else
					this_gap_len_b += 1
				end
			end
		end
		return scores
	end

	def trace_back(scores)
# FIXME cannot handle gap in last col or row
		aligned_seq_a = []
		aligned_seq_b = []
		alignment_score = 0

		if @penalize_endgaps then 
			i = @seq_a.size - 1
			j = @seq_b.size - 1
		else
			last_row = scores.last
			last_col = scores.collect {|row| row.last}
			i = last_col.each_with_index.max[1] # index of maximal value
			j = last_row.each_with_index.max[1]
		end

		this_gap_len_a, this_gap_len_b = 1, 1

		while ( i >= 0 || j >= 0) do 

			# require index really larger than 0
			# index = 0 will cause access of very last row/col ( due to rubies use of index (0-1) )

			if ( i > 0 && j > 0 && 
				scores[i][j] == scores[i-1][j-1] + aa_score(@seq_a[i], @seq_b[j]) ) then 

				aligned_seq_a.unshift(@seq_a[i])
				aligned_seq_b.unshift(@seq_b[j])

				alignment_score += aa_score(@seq_a[i], @seq_b[j])
				this_gap_len_a, this_gap_len_b = 1, 1

				i = i-1
				j = j-1

			elsif ( i > 0 && scores[i][j] == scores[i-1][j] + gap_cost(this_gap_len_b) )

				aligned_seq_a.unshift(@seq_a[i])
				aligned_seq_b.unshift("-")

				this_gap_len_b += 1
				alignment_score += gap_cost(this_gap_len_b)

				i = i-1

			elsif ( j > 0 && scores[i][j] == scores[i][j-1] + gap_cost(this_gap_len_a) )

				aligned_seq_a.unshift("-")
				aligned_seq_b.unshift(@seq_b[j])

				this_gap_len_a += 1
				alignment_score += gap_cost(this_gap_len_a)

				j = j-1

			elsif i * j == 0 
				# one or both values is zero
				# assemble remaining parts of both seqs, fill the other seq with gaps if neccessary

				max_n_missing_chars = [i,j].max

				aligned_seq_a.unshift( @seq_a[0..i].rjust(max_n_missing_chars, "-") )
				aligned_seq_b.unshift( @seq_b[0..j].rjust(max_n_missing_chars, "-") )

				if @penalize_endgaps then 
					alignment_score += gap_cost( max_n_missing_chars )
				end

				break

			else
				# this case should never happen
				Helper.abort "Needleman-Wunsch failed."
			end
					
		end

		max_score = scores.last.last

		return aligned_seq_a.join(""), aligned_seq_b.join(""), alignment_score.to_f/max_score
	end

	def index_of_max_2d_array(arr)
		max_val, ind_max_row, ind_max_col = -Float::INFINITY,0,0
		arr.each_with_index do |row, row_ind|
			this_max, col_ind = row.each_with_index.max
			if this_max > max_val then 
				max_val = this_max
				ind_max_row = row_ind
				ind_max_col = col_ind
			end
		end

		return ind_max_row, ind_max_col
	end

	def gap_cost(gap_len)
		if @lin_gap_pen then 
			# use linear gap penalty
			return lin_gap_pen(gap_len)
		else
			return affine_lin_gap_pen(gap_len)
		end
	end
	def lin_gap_pen(gap_len)
		gap_pen = 5 # cost of each gap
		return -(gap_len * gap_pen)
	end
	def affine_lin_gap_pen(gap_len)
		gap_open = 7 # cost to open a gap
		gap_ext = 1 # cost to extend a gap
		return -gap_open - (gap_len - 1) * gap_ext
	end

	def aa_score(a, b)
		@blosum62[(a.chr + b.chr).upcase].to_i
	end

	def read_in_blosum62
		matrix = {}

		file_path = File.join(File.dirname(__FILE__),'blosum62')
		Helper.file_exist_or_die(file_path)

		# convert blosum matrix to hash
		aas = nil # the order of amino acids in matrix

		IO.foreach(file_path) do |line|
			line.chomp!
			if line.start_with?("#") then 
				# comment line
				next
			elsif aas.nil?
				# definition of amino acid order
				aas = line.split(" ")
				next
			else
				# a matrix line
				cells = line.split(" ")
				aa1 = cells.delete_at(0)

				# convert matix to hash
				aas.each_with_index do |aa2, ind|
					matrix[aa1.to_s.upcase + aa2.to_s.upcase] = cells[ind]
				end
			end
		end

		return matrix
 	end

end