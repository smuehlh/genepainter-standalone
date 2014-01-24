module Helper
	extend self

	def open_log(file)
		@fh_log= File.open(file, "w")
		log "Logfile created on #{Time.now} by GenePainter"
	end
	def log(msg)
		@fh_log.puts msg.to_s
	end
	def close_log
		@fh_log.close
	end

	def abort(msg)
		$stderr.puts "\nFatal error: #{msg.to_s}"
		exit 1
	end

	def warn(msg)
		$stderr.puts msg.to_s
	end

	def dir_exist_or_die(path)
		if ! FileTest.directory?(path) then
			abort "Directory #{path} does not exist."
		end
	end

	def file_exist_or_die(path)
		if ! FileTest.file?(path) then
			abort "File #{path} does not exist."
		end
	end

	# inform which data are not used
	def print_intersect_and_diff_between_alignment_and_gene(alignment_names, gene_names)
		common = alignment_names & gene_names
		missing_genes = alignment_names - gene_names
		missing_seqs = gene_names - alignment_names
		puts ""
		if missing_genes.any? then
			puts "No gene structures for: #{missing_genes.join(", ")}. Ignore them.\n"
		end
		if missing_seqs.any? then
			puts "No aligned sequences for: #{missing_seqs.join(", ")}. Ignore them.\n"
		end
		puts "Using #{common.size} sequences."

		# write to log file which seqs and genes exactly are used
		log ""
		log "Used #{common.join(", ")} for computation."

	end

	# sub_str can be either a string or a regex
	def find_each_index(search_str,sub_str)
		i = -1
		inds = []
		while i = search_str.index(sub_str, i+1)
			inds << i
		end
		return inds
	end

	def convert_number_to_human_readable_string(num)
		# i.e. add thousand separator
		num = num.to_s 
		# 1) split the string in chunks of 3 [scan] - but start from last position [reverse]
		# 2) join the chunks by a comma [join] - but reverse again to see orignal number [reverse]
		return num.reverse.scan(/.{1,3}/).join(",").reverse
	end

	def ruby2human_counting(num)
		num + 1
	end
end

class Array
	# multiplies each element in array with x
	# recursively, if array elements are arrays itself
	def multiply_by(x)
		collect do |v|
			case(v)
			when Array
				# if item in array is an array itself, apply same method on it
				v.multiply_by(x)
			else
				v * x
			end
		end
	end

	def sum
		self.inject(0, :+)
	end
end
