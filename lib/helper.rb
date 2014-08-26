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

	def log_error(function_error_was_raised, msg)
		$stderr.puts "\nError while executing '#{function_error_was_raised}'"
		log "Error while executing '#{function_error_was_raised}': #{msg}"
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

	def dir_writable_or_die(path)
		if ! File.writable?(File.expand_path(path)) then 
			abort "Directory #{path} is not writable."
		end
	end

	def file_exist_or_die(path)
		if ! FileTest.file?(path) then
			abort "File #{path} does not exist."
		end
	end

	def file_exist?(path)
		FileTest.file?(path)
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
		log "Using #{common.join(", ")} genes for computation."

	end

	# inform which data are used and not used for taxonomy
	def print_intersect_and_diff_between_taxonomy_genes_and_selected_taxa(all_genes, genes_with_taxonomy, selected_taxa, genes_within_selected_taxa)
		missing_tax = all_genes - genes_with_taxonomy
		puts ""

		# simultaniously write to log 
		log ""
		if missing_tax.any? then
			puts "No taxonomy for: #{missing_tax.join(", ")}."
		end
		if genes_with_taxonomy.any? then 
			puts "Taxonomy for: #{genes_with_taxonomy.join(", ")}."
			log "Genes with taxonomic information: #{genes_with_taxonomy.join(", ")}."
		end
		if selected_taxa.any? then 
			if genes_within_selected_taxa.any? then
				puts "Genes belonging to selected taxa >#{selected_taxa.join(", ")}<: #{genes_within_selected_taxa.join(", ")}"
				log "Genes belonging to selected taxa #{selected_taxa.join(", ")}: #{genes_within_selected_taxa.join(", ")}"
			else
				puts "No genes belonging to selected taxa >#{selected_taxa.join(", ")}<"
				log "No genes belonging to selected taxa #{selected_taxa.join(", ")}"
			end
		end
	end

	# inform which data are used after selection is applied 
	def print_selection(gene_names, selection_type, selection_criterium)
		puts ""
		puts "Appling selection to #{selection_type} based on #{selection_criterium}."
		puts "Restricting data set to #{gene_names.join(", ")}."
		puts "Using #{gene_names.size} sequences."

		log ""
		log "Appling selection to #{selection_type} based on #{selection_criterium}."
		log "Using #{gene_names.join(", ")} genes for computation."
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
	def human2ruby_counting(num)
		num - 1
	end

	def word_frequency(arr)
		res = Hash.new(0)
		arr.each { |a| res[a] += 1 }
		res.delete(nil)
		return res
	end

	def sanitize_taxon_name(str)
		str = str.gsub(/\s/, '.')
		return str.gsub('[^A-Za-z0-9_\.\-]', '').gsub('..', '.')
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

	# mimics set method "subset?"
	def is_subset?(other_arr)
		return false if other_arr.size < size
		all? { |obj| other_arr.include?(obj) }
	end
	# mimics set method "disjoint?"
	def is_disjoint_set?(other_arr)
		! is_overlapping_set?(other_arr)
	end
	# mimics set method "intersect?"
	def is_overlapping_set?(other_arr)
		# (self & other_arr).any?
		if size < other_arr.size then 
			any? { |obj| other_arr.include?(obj) }
		else
			other_arr.any? { |obj| include?(obj) }
		end
	end
	def intersection(other_arr)
		self & other_arr
	end
end
