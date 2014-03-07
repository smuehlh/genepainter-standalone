class TaxonomyToGene

	attr_reader :known_genes, :taxa_with_tax_obj, :species_with_corresponding_genes
	
	# path_to_taxdump: path to NCBI taxonomy dump file
	# path_to_linked_list: path to file specifing which fasta header belongs to which species
	def initialize(path_to_taxdump, path_to_linked_list_genenames_speciesnames)

		@grep = `which grep`.chomp

		@species_with_corresponding_genes = map_genenames_to_speciesnames( path_to_linked_list_genenames_speciesnames )
		@known_genes = @species_with_corresponding_genes.values.flatten.uniq

		@taxa_with_tax_obj = extract_taxonomy_of_species_in_alignment_from_taxdump( path_to_taxdump ) # Hash: key= taxon, value: taxonomie-object
	end

	# read in file containing mapping between fasta header (=gene names) and species
	def map_genenames_to_speciesnames(path_to_linked_list)
		species2genes = {}
		IO.foreach(path_to_linked_list) do |line|
			line = line.chomp
			genes, species = line.split(/:/x) # ignore white spaces surrounding ":"
			if genes.nil? || species.nil? then
				Helper.abort "Invalid syntax in file #{path_to_linked_list}. Expecting \':\'-separated list of genes and species"
			end
			species = species.strip # remove leading & trailing white spaces
			genes = genes.split(/,/).map { |g| g.strip }
			species2genes[species] = genes
		end
		return species2genes
	end

	# extract lineage for all species present in alignment from taxdump
	# uses grep, as taxdump is large file and ram is restricted
	def extract_taxonomy_of_species_in_alignment_from_taxdump( path_to_taxdump )

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
		threads = []
		get_all_species_linked_to_genes.each do |species|

			threads << Thread.new{ Thread.current[:output] = get_lineage_by_species(species, path_to_nodesdmp, path_to_namesdmp) }

		end
		threads.each do |tr|
			tr.join
			species, lineage = tr[:output]
			lineage_by_species[species] = lineage
		end

		# convert lineages to taxonomy objects (cant do this in one step, as distance to root is only known when having complete lingeage)
		tax_objs_by_name = {}

		get_all_species_linked_to_genes.combination(2).each do |species1, species2|

			lineage1 = lineage_by_species[species1].reverse
			lineage2 = lineage_by_species[species2].reverse

			lca = (lineage1 & lineage2).last # last as lineages are in reverse order: from root to species

			add_or_update_taxonomy_obj(tax_objs_by_name, lineage1, species1, lca)
			add_or_update_taxonomy_obj(tax_objs_by_name, lineage2, species2, lca)

		end

		return tax_objs_by_name
	end

	# creates/updates taxonomy object for every taxon of lineage
	# input lineage is in reverse order: from root to species
	def add_or_update_taxonomy_obj(tax_objs_by_name, lineage, species, lca)

		lineage.each_with_index do |taxon, ind|
			if tax_objs_by_name[taxon] then 
				# update children
				tax_objs_by_name[taxon].add_child(species)
				
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
					descendant, # child (node)
					ind, # distance to root
					species # child (leaf)
				)
			end
			if taxon == lca then 
				tax_objs_by_name[taxon].add_lca(species)
			end
		end
	end

	# returns to which species (list) the input genes belong
	def get_species_by_genes(genes_list)
		species_list = []
		@species_with_corresponding_genes.each do |species, genes|
			if genes.is_overlapping_set?(genes_list) then 
				species_list.push( species )
			end
		end
		return species_list
	end

	# finds the last common ancestor of species defined species_list in putative_last_common_ancestor_list
	def get_last_common_ancestor_of(species_list)
		
		lca = ""
		@taxa_with_tax_obj.each do |taxon, taxon_obj|
			# taxon is lca of all species in list
			if taxon_obj.is_last_common_ancestor? && species_list.is_subset?( @taxa_with_tax_obj[taxon].last_common_ancestor_of ) then 
				if lca.empty? then 
					# this is the first lca found, save it
					lca = taxon
				else
					# compare this lca to last found lca: chose lca with greater distance to root
					distance_to_root_this_lca = @taxa_with_tax_obj[taxon].distance_to_root
					distance_to_root_last_found_lca = @taxa_with_tax_obj[lca].distance_to_root
					if distance_to_root_this_lca > distance_to_root_last_found_lca then 
						lca = taxon
					end
				end
			end
		end
		return lca
	end

	# finds first uniq ancestor for every species associated with gene from genes list
	# also returns number of genes per first uniq ancestor
	def get_first_uniq_ancestors_with_frequencies_by_genes(genes_list, lca)

		first_uniq_list = []
		occurence_first_uniq = []

		children_of_lca = @taxa_with_tax_obj[lca].children
		descendants_of_lca = @taxa_with_tax_obj[lca].descendants

		@species_with_corresponding_genes.each do |species, genes|
			if genes.is_overlapping_set?(genes_list) then

				n_genes_encoded_by_this_species = genes.intersection(genes_list).size

				# first uniq ancestor of this species
				first_uniq = descendants_of_lca.find do |taxon|
					@taxa_with_tax_obj[taxon].children.include?(species)
				end
				if first_uniq then 

					# if lca is a species, then there is no first uniq ancestor
					ind_in_results = first_uniq_list.index(first_uniq)
					if ind_in_results then
						# update counts
						occurence_first_uniq[ind_in_results] += n_genes_encoded_by_this_species
					else
						# add first_uniq and update counts
						first_uniq_list.push(first_uniq)
						occurence_first_uniq.push(n_genes_encoded_by_this_species)
					end
				end
			end
		end 

		return first_uniq_list, occurence_first_uniq

	end

	# # finds first uniq ancestors of species list
	# # input: species list
	# # species frequency (same order as species list)
	# # if optional argument lca ist not provided, the last common ancestor of all species will be calculated first
	# # output might be empty, if lca is a species ( a species has no further descendants)
	# def get_first_uniq_ancestors_with_frequency_of(species_list, species_frequencies, *lca)
	# 	first_uniq_list = []
	# 	occurence_first_uniq = []
	# 	if lca.empty? then 
	# 		lca = get_last_common_ancestor_of(species_list)
	# 	else
	# 		lca = lca.first
	# 	end
	# 	@taxa_with_tax_obj[lca].descendants.each do |taxon|
	# 		children = @taxa_with_tax_obj[taxon].children
	# 		if children.is_overlapping_set?(species_list) then 
	# 			first_uniq_list.push( taxon )
	# 			occurence_first_uniq.push( children.intersection(species_list).size )
	# 		end
	# 	end
	# 	return first_uniq_list, occurence_first_uniq
	# end

	# returns genes encoded by children of input
	def get_genes_encoded_by_taxon(taxon)
		if @taxa_with_tax_obj[taxon] then 
			return @taxa_with_tax_obj[taxon].children.collect do |child|
				get_genes_encoded_by_species(child)
			end
		else
			return []
		end
	end
	def get_genes_encoded_by_species(species)
		if @species_with_corresponding_genes[species] then 
			return @species_with_corresponding_genes[species].uniq
		else
			return []
		end
	end

	def get_all_species_linked_to_genes
		@species_with_corresponding_genes.keys
	end

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
		output = io.read # returns string,buffer or nil
		io.close
		if output.respond_to?('split') then 
			output.split(/\t?\|\t?/)
		else
			return ["",""] 
		end
	end
end
