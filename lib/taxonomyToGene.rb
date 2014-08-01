class TaxonomyToGene

	attr_reader :known_genes, :taxa_with_tax_obj, :species_with_corresponding_genes, 
		:species_not_found_in_NCBI_tax, :species_without_gene_structures

	# path_to_taxdump: path to NCBI taxonomy dump file
	# path_to_linked_list: path to file specifing which fasta header belongs to which species
	# genes_with_data: list of all genes with gene structure (ignore all other genes named in path_to_linked_list)
	def initialize(path_to_taxonomy, path_to_linked_list_genenames_speciesnames, is_no_grep, is_call_grep_with_nice, genes_with_data )

		@species_not_found_in_NCBI_tax = []
		@species_without_gene_structures = []

		@grep = `which grep`.chomp
		if is_call_grep_with_nice then 
			@be_nice = true
		end

		@species_with_corresponding_genes = self.class.map_genenames_to_speciesnames( path_to_linked_list_genenames_speciesnames, genes_with_data )
		@known_genes = @species_with_corresponding_genes.values.flatten.uniq
		
		species_with_corresponding_lineages = extract_taxonomy_for_species(path_to_taxonomy, is_no_grep)
		if @known_genes.empty? || species_with_corresponding_lineages.empty? then 
			Helper.abort "Mapping of genes to corresponding taxonomic lineage failed. No gene could be mapped."
		end

		@taxa_with_tax_obj = convert_lineages_to_taxonomy_objs(species_with_corresponding_lineages) # Hash: key= taxon, value: taxonomie-object
	end

	# include taxonomic info in gene obj: full lineage of corresponding species and index of the first uniq ancestor of that species
	def to_gene(gene_name)
		tax_info_for_gene = {lineage: [], ind_first_uniq: 0}
		if @known_genes.include?(gene_name) then 
			# connection gene-species is known
			species = get_species_by_gene(gene_name)
			lineage, ind_first_uniq_taxon = extract_full_lineage_and_first_uniq_ancestor_of(species)

			tax_info_for_gene[:lineage] = lineage
			tax_info_for_gene[:ind_first_uniq] = ind_first_uniq_taxon
		else
			# no connection to species present for this gene
		end
		return tax_info_for_gene
	end

	# returns to which species the gene belongs
	def get_species_by_gene(gene)
		@species_with_corresponding_genes.each do |species, genes_of_species|
			if genes_of_species.include?(gene) then 
				return species
			end
		end
		return ""
	end
	def extract_full_lineage_and_first_uniq_ancestor_of(species)
		lineage = []
		ind = 0
		if ! @taxa_with_tax_obj.has_key?(species) then 
			return lineage, ind
		end
		tax_obj = @taxa_with_tax_obj[species]
		species_name = tax_obj.name
		first_uniq_ancestor_name = nil
		while ! tax_obj.ancestor.empty? do

			# collect lineage
			lineage.unshift(tax_obj.name)
			tax_obj = @taxa_with_tax_obj[tax_obj.ancestor]

			# search for the first uniq taxon of this species
			if tax_obj.is_last_common_ancestor_of?(species) && first_uniq_ancestor_name.nil? then 
				# the descendant of the LCA was just added to the lineage...
				first_uniq_ancestor_name = lineage[0]
			end
		end
		lineage.unshift(tax_obj.name)
		ind = lineage.index(first_uniq_ancestor_name)
		return lineage, ind
	end

	# issue an error while parsing user-provided species or taxonomy file
	def self.error_while_parsing_file(file, line, msg)
		Helper.abort "Invalid syntax in file #{file} in line #{line}. #{msg}"
	end

	# read in file containing mapping between fasta header (=gene names) and species
	def self.map_genenames_to_speciesnames(path, mandatory_genes)
		species_with_genes = {}

		IO.foreach(path) do |line|
			line = line.chomp
			next if line.empty? # skip empty lines

			match_data = line.match(
				/^ 
				([^\"]+) # gene list
				\s*:\s* # delimiter between genes and species, might be surrounded by white spaces
				\"([^\"]+)\" # species, must be enclosed by double quotes
				$ # the end of line, to ensure that species does not contain any double quotes in name
				/x) # x: allows to comment on individual parts of regex

			if match_data.nil? || match_data.size != 3 then 
				error_while_parsing_file(path, line, 
					"Expected colon-separated list of gene(s) and species. Species must be enclosed in double quotes (\")"
				)
			end
			genes = match_data[1].split(/[;,]/).map { |g| g.strip } # remove leading and trailing white spaces
			species = match_data[2].strip.capitalize # remove leading and trailing white spaces, capitalize species name

			if ! genes.is_overlapping_set?(mandatory_genes) then 
				# for this genes, no gene structures exist. do not bother with their taxonomy
				@species_without_gene_structures |= [species]
				next
			end

			# collect mapping between species and genes
			if ! species_with_genes[species] then 
				species_with_genes[species] = []
			end
			species_with_genes[species] |= genes

		end

		return species_with_genes
	end

	def extract_taxonomy_for_species(path_to_tax, is_no_grep)
		if check_if_taxonomy_is_taxdump_or_list_of_lineages(path_to_tax) then 
			# tax file is NCBI taxonomy dump
			return extract_taxonomy_of_species_in_alignment_from_taxdump( path_to_tax, is_no_grep ) 
		else
			# tax file is no taxonomy dump, but contains only lineage
			return extract_taxonomy_of_species_in_alignment_from_taxonomylist(path_to_tax)
		end
	end
	def check_if_taxonomy_is_taxdump_or_list_of_lineages(path)
		dir_taxdump = File.dirname(path)
		tar = `which tar`.chomp
		is_success = system("cd #{dir_taxdump} && #{tar} -tf #{File.absolute_path(path)} >/dev/null 2>&1")
		return is_success
	end

	# extract lineage for all species present in alignment from taxdump
	# uses grep, as taxdump is large file and ram is restricted
	def extract_taxonomy_of_species_in_alignment_from_taxdump( path_to_taxdump, is_no_grep )

		lineage_by_species = {}

		# extract names and nodes from tar archive
		dir_taxdump = File.dirname(path_to_taxdump)
		path_to_namesdmp = File.join( dir_taxdump, "names.dmp" )
		path_to_nodesdmp = File.join( dir_taxdump, "nodes.dmp")
		if ! ( Helper.file_exist?( path_to_namesdmp ) && Helper.file_exist?( path_to_nodesdmp ) ) then 
			# last attemp to get files; maybe they simply need to be extracted
			tar = `which tar`.chomp
			is_success = system("cd #{dir_taxdump} && #{tar} -xvf #{File.absolute_path(path_to_taxdump)} names.dmp nodes.dmp >/dev/null 2>&1")
			Helper.file_exist_or_die(path_to_nodesdmp)
			Helper.file_exist_or_die(path_to_namesdmp)

		end

		# process all species names and find their taxid

		if is_no_grep then 
			# read names and nodes into RAM

			nodes_with_parents_and_name = {}
			species_with_taxid = {}

			IO.foreach(path_to_nodesdmp) do |line|
				line.chomp!
				parts = line.split("\t|\t")
				taxid = parts[0].to_i
				parent_taxid = parts[1].to_i
				nodes_with_parents_and_name[taxid] = [ parent_taxid ]
			end

			IO.foreach(path_to_namesdmp) do |line|
				line.chomp!
				if line.include?("scientific name") then 
					parts = line.split("\t|\t")
					taxid = parts[0].to_i
					name = parts[1]
					nodes_with_parents_and_name[taxid].push(name)
					if get_all_species_linked_to_genes.include?(name) then 
						species_with_taxid[name] = taxid
					end
				end
			end

			get_all_species_linked_to_genes.each do |species|

				species, lineage = get_lineage_by_species_no_grep(species, nodes_with_parents_and_name, species_with_taxid)
				# parsing taxonomy for species might have been not successful!
				# in that case, nil values will be returned
				if lineage then 
					lineage_by_species[species] = lineage
				else
					# remove species from list of linked species
					# collect species name to output information
					unlink_species_from_gene(species)
					@species_not_found_in_NCBI_tax.push( species )
				end	
			end

		else
			get_all_species_linked_to_genes.each_slice(10) do |slice|
				print "."

				threads = []
				slice.each do |species|
					threads << Thread.new{ Thread.current[:output] = get_lineage_by_species(species, path_to_nodesdmp, path_to_namesdmp) }
				end
				threads.each do |tr|
					tr.join
					species, lineage = tr[:output]
					# parsing taxonomy for species might have been not successful!
					# in that case, nil values will be returned
					if lineage then 
						lineage_by_species[species] = lineage
					else
						# remove species from list of linked species
						# collect species name to output information
						unlink_species_from_gene(species)
						@species_not_found_in_NCBI_tax.push( species )
					end
				end
			end

		end

		return lineage_by_species
	end

	def extract_taxonomy_of_species_in_alignment_from_taxonomylist(path)
		lineage_by_species = {}

		species_with_lineage_without_genes = []

		IO.foreach(path) do |line|
			line = line.chomp
			next if line.empty? # skip empty lines

			lineage = line.split(";").map { |g| g.strip.capitalize } # remove leading and trailing white spaces, capitalize all taxa
			if lineage.size < 2 then 
				# lineage must consist at least of root and species
				self.class.error_while_parsing_file(path, line, "Expected semicolon-separated list of taxa from root to species.")
			end
			species = lineage[-1] # last taxon listed must be species
			if @species_with_corresponding_genes[species] then 
				# genes are mapped to this species, so collect the taxonomy

				lineage_by_species[species] = lineage.reverse
			else
				species_with_lineage_without_genes |= [species]
			end
		end

		# handle species with taxonomic information, but without genes
		if species_with_lineage_without_genes.any? then 
			Helper.warn "Found taxonomic lineage for species that were not mapped to any genes: #{species_with_lineage_without_genes.join(", ")}. Ignore them."
		end

		# handle species without taxonomic information, but with genes
		(get_all_species_linked_to_genes - lineage_by_species.keys).each do |species|
			# remove species from list of linked species
			# collect species name to output information
			unlink_species_from_gene(species)
			@species_not_found_in_NCBI_tax.push( species )
		end
		return lineage_by_species
	end

	def convert_lineages_to_taxonomy_objs(lineage_by_species)
		# convert lineages to taxonomy objects (cant do this in one step, as distance to root is only known when having complete lingeage)
		tax_objs_by_name = {}

		get_all_species_linked_to_genes.combination(2).each do |species1, species2|

			lineage1 = lineage_by_species[species1].reverse
			lineage2 = lineage_by_species[species2].reverse

			lca = (lineage1 & lineage2).last # last as lineages are in reverse order: from root to species

			self.class.add_or_update_taxonomy_obj(tax_objs_by_name, lineage1, species1, lca)
			self.class.add_or_update_taxonomy_obj(tax_objs_by_name, lineage2, species2, lca)

		end

		return tax_objs_by_name
	end

	# creates/updates taxonomy object for every taxon of lineage
	# input lineage is in reverse order: from root to species
	def self.add_or_update_taxonomy_obj(tax_objs_by_name, lineage, species, lca)

		lineage.each_with_index do |taxon, ind|
			if tax_objs_by_name[taxon] then 
			
				# update descendants
				if ind != (lineage.size - 1 ) then 
					descendant = lineage[ind+1]
					tax_objs_by_name[taxon].add_descendant(descendant)
				end
			else
				# create object

				# find ancestor and descendant of taxon
				if ind == 0 then 
					ancestor = ""
				else
					ancestor = lineage[ind-1]
				end
				if ind == (lineage.size - 1) then 
					descendant = ""
				else
					descendant = lineage[ind+1]
				end

				tax_objs_by_name[taxon] = Taxonomy.new(
					taxon, # name
					ancestor,
					descendant # child (node)
				)
			end
			if taxon == lca then 
				tax_objs_by_name[taxon].add_lca(species)
			end
		end
	end

	def get_all_species_linked_to_genes
		@species_with_corresponding_genes.keys
	end

	def unlink_species_from_gene(species)
		# delete species and its corresponding genes
		genes_belonging_to_species = @species_with_corresponding_genes.delete(species) # delete returns what is being deleted
		@known_genes = @known_genes - genes_belonging_to_species
	end

	### methods to extract lineage via grep ###
	def get_lineage_by_species(species, path_to_nodesdmp, path_to_namesdmp)
		# get taxonomy of this species
		lineage = [species]

		taxid = get_taxid_by_name(species, path_to_namesdmp)
		while taxid != "1" do
			parent_taxid, parent_name = get_parent_taxid_with_name(taxid, path_to_nodesdmp, path_to_namesdmp)
			lineage << parent_name.capitalize
			taxid = parent_taxid
		end

		return species, lineage
	rescue
		return species, nil
	end

	def get_lineage_by_species_no_grep(species, nodes_with_parents_and_name, species_with_taxid)
		# get taxonomy of this species
		lineage = [species]

		taxid = species_with_taxid[species]
		while taxid != 1 do
			parent_taxid, parent_name = nodes_with_parents_and_name[taxid]
			lineage << parent_name.capitalize
			taxid = parent_taxid
		end

		return species, lineage
	rescue
		return species, nil
	end

	def get_parent_taxid_with_name(taxid, path_to_nodesdmp, path_to_namesdmp)
		# find parent tax id
		parent_taxid = get_parent_taxid(taxid, path_to_nodesdmp)

		# find corresponding name
		parent_name = get_name_by_taxid(parent_taxid, path_to_namesdmp)

		return parent_taxid, parent_name
	end

	# structure of names.dmp
	# tax_id\t|\tname\t|\t ...
	# structure of nodes.dmp
	# tax_id\t|\tparent tax_id\t|\t ...
	def get_parent_taxid(taxid, file)
		output_parts = grep_in_dmp("^#{taxid}\t", file)
		return output_parts[1]
	end
	def get_name_by_taxid(taxid, file)
		# output_parts = grep_in_dmp("^#{taxid}\t", file)
		output_parts = grep_pipe_without_maxcount("^#{taxid}\t", "scientific name", file)
		return output_parts[1]
	end
	def get_taxid_by_name(name, file)
		output_parts = grep_in_dmp("\t#{name}\t", file)
		return output_parts[0]
	end
	def grep_in_dmp(search_term, file)
		# consider wrapping search_term in "\t" to avoid false positives
		command = "#{@grep} --max-count=1 -i -P '#{search_term}' #{file}"
		if @be_nice then 
			command = "nice #{command}"
		end
		io = IO.popen(command)
		output = io.read # returns string,buffer or nil
		io.close
		if output.respond_to?('split') then 
			output.split(/\t?\|\t?/)
		else
			return ["",""] 
		end
	end
	def grep_pipe_without_maxcount(search_term1, search_term2, file)
		command = "#{@grep} -i -P '#{search_term1}' #{file} | #{@grep} -P '#{search_term2}'"
		if @be_nice then 
			command = "nice #{command}"
		end
		io = IO.popen(command)
		output = io.read # returns string,buffer or nil
		io.close
		if output.respond_to?('split') then 
			output.split(/\t?\|\t?/)
		else
			return ["",""] 
		end
	end

end
