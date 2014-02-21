class Taxonomy

	# path_to_taxdump: path to NCBI taxonomy dump file
	# path_to_linked_list: path to file specifing which fasta header belongs to which species
	def initialize(path_to_taxdump, path_to_linked_list)

		@grep = `which grep`.chomp

		@gene_names2species_names = map_genenames_to_species(path_to_linked_list) # key: gene (=fasta header; gene.name), value: species name
		@species_with_parent_taxa = reduce_taxdump_to_taxa_from_alignment(path_to_taxdump) # key: species, value: array of ranks upto root of taxonomy 
	end

	def reduce_taxdump_to_taxa_from_alignment(path_to_taxdump)
		# due to RAM: parse NCBI tax dump only with grep; which can handle gzipped-file

		taxa_with_parents = {}

		# extract names and nodes from tar archive
		dir_taxdump = File.dirname(path_to_taxdump)
		path_to_namesdmp = File.join( dir_taxdump, "names.dmp" )
		path_to_nodesdmp = File.join( dir_taxdump, "nodes.dmp")
		if ! ( Helper.file_exist?( path_to_namesdmp ) && Helper.file_exist?( path_to_nodesdmp ) ) then 
			tar = `which tar`.chomp
			is_success = system("cd #{dir_taxdump} && #{tar} -xvf #{File.absolute_path(path_to_taxdump)} names.dmp nodes.dmp >/dev/null 2>&1")
			Helper.file_exist_or_die(path_to_nodesdmp)
			Helper.file_exist_or_die(path_to_namesdmp)
		end

		# process all species names and find their taxid
		@gene_names2species_names.values.uniq.each do |species|
			taxid = get_taxid_by_name(species, path_to_namesdmp)

			if ! taxa_with_parents[species] then 
				# never obtained taxonomy for this species before
				taxa_with_parents[species] = []
				while taxid != "1" do 
					parent_taxid, parent_name = get_parent_taxid_with_name(taxid, path_to_nodesdmp, path_to_namesdmp)
					taxa_with_parents[species] << parent_name.capitalize
					taxid = parent_taxid
				end
			end
		end
		return taxa_with_parents
	end

	def map_genenames_to_species(path_to_linked_list)
		genes2names = {}
		IO.foreach(path_to_linked_list) do |line|
			line = line.chomp
			genes, species = line.split(/:/x) # ignore white spaces surrounding ":"
			if genes.nil? || species.nil? then
				Helper.abort "Invalid syntax in file #{path_to_linked_list}. Expecting \':\'-separated list of genes and species"
			end
			species = species.strip # remove leading & trailing white spaces
			genes.split(/,/).each do |gene|
				# genes might be ","-separated list
				gene = gene.strip 
				genes2names[gene] = species
			end
		end
		return genes2names
	end

	# make sure in pre-processing, that both gene_names and selected_taxa are all capitalized
	# otherwise, might miss some matches due to uppercase-lowercase problems
	def get_genenames_belonging_to_selected_taxa(gene_names, selected_taxa)

		reduced_names_list = []

		# iterate through all genes for which 
		gene_names.each do |gene|
			if species = @gene_names2species_names[gene] then 
				parents = @species_with_parent_taxa[species]
				if ( parents & selected_taxa ).any? then 
					reduced_names_list << gene
				else
					# no match, skip
				end
			else
				# no species entry
				# skip
			end
		end
		return reduced_names_list
	end

	def get_all_parents_for_genes(gene_names)
		parents = {} # array of parents for each gene
		gene_names.each do |gene|
			parents[gene] = []

			if species = @gene_names2species_names[gene] then 
				parents[gene] = @species_with_parent_taxa[species]
			else
				# no species entry
				# skip
			end
		end
		return parents
	end

	def get_all_inner_nodes_with_children(gene_names)
		nodes = {} # array of children (= genes) for each parent

		# TODO vergiss all die inner nodes, wo sich "nichts tut"; mit first_uniq?

		gene_names.combination(2).each do |gene1, gene2|

			if species1 = @gene_names2species_names[gene1] && species2 = @gene_names2species_names[gene2] then 
				parents1 = @species_with_parent_taxa[species1]
				parents2 = @species_with_parent_taxa[species2]

				first_uniq = self.class.first_uniq_ancestors_of_genes( [parents1, parents2] )
puts first_uniq
				debugger
				puts "hm,hm"
			else
				# skip
			end 
		end
	end

	# return the last common ancestor
	# expects array of arrays, order: species to root
	# returns string with last common ancestor or empty string
	def self.last_common_ancestor_of_genes(ancestors_of_genes)

		# build intersect of all ancestors
		# first element of list is last common ancestor (the most recent of all common ancestors)
		all_common_ancestors = ancestors_of_genes.inject{ |ancestors, intersect| ancestors & intersect }
		return all_common_ancestors.first 
	rescue NoMethodError
		return ""
	end

	# returns first ancestor of every gene which follows the last common ancestor of all genes
	# expects array of arrays, order: species to root
	# returns array of strings, will be smaller than array of genes if lca is last ancestor at all (e.g. before species-level)
	def self.first_uniq_ancestors_of_genes(ancestors_of_genes)

		lca = last_common_ancestor_of_genes(ancestors_of_genes)

		first_uniq_ancestor_list = []

		ancestors_of_genes.each do |ancestors|
			ind = ancestors.index(lca)

			# convert index of last common to index of first uniq ancestor
			if ind > 0 then 
				ind_first_after_lca = ind - 1 # -1 because list goes from recent to most ancient
				first_uniq_ancestor_list.push( ancestors[ind_first_after_lca] )
			end
		end
		return first_uniq_ancestor_list
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
		output_parts = grep_in_dmp("^#{taxid}\t", file)
		return output_parts[1]
	end
	def get_taxid_by_name(name, file)
		output_parts = grep_in_dmp("\t#{name}\t", file)
		return output_parts[0]
	end
	def grep_in_dmp(search_term, file)
		# consider wrapping search_term in "\t" to avoid false positives
		io = IO.popen([@grep, "--max-count=1", "-i", "-P", "#{search_term}", file])
		output = io.readlines
		io.close
		return output.first.split(/\t?\|\t?/)
	end
end
