module Helper
	extend self

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
	def print_intersect_and_diff_between_alignment_and_gene(alignment_names, gene_names, fh_log)
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
		fh_log.puts ""
		fh_log.puts "Used #{common.join(", ")} for computation."

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
end
