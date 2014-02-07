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
