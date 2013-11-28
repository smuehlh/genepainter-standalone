module Helper
	extend self

	def abort(msg)
		$stderr.puts "\nFatal error: #{msg.to_s}"
		exit 1
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

end