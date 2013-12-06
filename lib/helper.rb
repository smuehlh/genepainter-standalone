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

	module SvgPainter
		extend self
		def line(x1, y1, x2, y2, color)
			return "<line x1=\"#{x1}\" "\
				"y1=\"#{y1}\" "\
				"x2=\"#{x2}\" "\
				"y2=\"#{y2}\" "\
				"stroke=\"#{color}\" "\
				"/>"
		end

		def box(x, y, width, height, color)
			return "<rect x=\"#{x}\" "\
				"y=\"#{y}\" "\
				"width=\"#{width}\" "\
				"height=\"#{height}\" "\
				"fill=\"#{color}\" "\
				"/>"
		end

		def text(x,y,txt)
			return "<text x=\"#{x}\" "\
				"y=\"#{y}\"> "\
				"#{txt} "\
				"</text>"
		end

		def header(width=1000,height=500)
			return "<svg version=\"1.1\" "\
				"baseProfile=\"full\" "\
				"width=\"#{width}\" "\
				"height=\"#{height}\" "\
				"viewBox=\"0 0 #{width.to_f/10} #{height.to_f/10}\" "\
				"xmlns=\"http://www.w3.org/2000/svg\">"
		end
		def footer
			return "</svg>"
		end
	end

end
