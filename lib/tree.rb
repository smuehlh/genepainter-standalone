class Tree

	def initialize(all_lineages, last_common_ancestor, mandatory_taxa)
		# create all taxonomy objects needed for tree
		build_tree(all_lineages, last_common_ancestor, mandatory_taxa) 
	end
	def build_tree(all_lineages, last_common_ancestor, mandatory_taxa)
		@taxa_with_tax_objs = {}

		add_root_to_tree(last_common_ancestor)

		add_mandatory_nodes_to_tree(all_lineages, mandatory_taxa)

		add_missing_species_to_tree_fewest_needed_leaves(all_lineages)
		# add_missing_species_to_tree(all_lineages, first_uniq_ancestors_by_species)

		validate_tree(all_lineages)

		annotate_nodes_with_leaves_below
	end

	def get_all_nodes
		@taxa_with_tax_objs.keys
	end

	def is_leaf?(node)
		@taxa_with_tax_objs[node].is_leaf?
	end

	def get_leaves
		return @taxa_with_tax_objs.collect{|k,v| k if v.is_leaf?}.compact
	end

	# trust method that called initialize() that the root indeed occurs in all lineages
	# just in case that nothing occurs in all lineage (which should never happen), invent a root
	def add_root_to_tree(root)
		if ! root then
			Helper.warn "No common ancestor of all lineages found when building tree."
			root = "root"
		end
		@taxa_with_tax_objs[root] = create_dummy_taxonomy_obj(root)
		@root = root
	end
	def add_mandatory_nodes_to_tree(all_lineages, mandatory_taxa)

		mandatory_taxa.each do |current_node|

			if @taxa_with_tax_objs[current_node] then 
				# current node is already part of the tree
				next
			end
			lineage_root_to_current_node = extract_lineage_root_to_node(current_node, all_lineages)

			parent_node = find_parent_in_tree(current_node, lineage_root_to_current_node)

			if is_no_cycle_between_parent_and_node(parent_node, current_node) then 
				@taxa_with_tax_objs[current_node] = create_dummy_taxonomy_obj(current_node)
				link_node_to_parent(current_node, parent_node)
			end
			verify_links_siblings_to_parent(parent_node, current_node, lineage_root_to_current_node, all_lineages)

		end # each mandatory taxon
	end

	def add_missing_species_to_tree(all_lineages, first_uniq_ancestors_by_species)
		all_lineages.each do |lineage_root_to_species|
			current_species = lineage_root_to_species.last

			parent_node = find_parent_in_tree(current_species, lineage_root_to_species)

			if ! @taxa_with_tax_objs[parent_node].is_leaf? then 
				current_node = first_uniq_ancestors_by_species[current_species]

				if is_no_cycle_between_parent_and_node(parent_node, current_node) then 
					@taxa_with_tax_objs[current_node] = create_dummy_taxonomy_obj(current_node)
					link_node_to_parent(current_node, parent_node)
				end
				verify_links_siblings_to_parent(parent_node, current_node, lineage_root_to_species, all_lineages)

			end # if parent_node is no leaf

		end # each species
	end

	def add_missing_species_to_tree_fewest_needed_leaves(all_lineages)
		all_lineages.each do |lineage_root_to_species|
			current_species = lineage_root_to_species.last

			# is species itself part of the tree? (and hopefully a leaf)?
			if @taxa_with_tax_objs[current_species] && @taxa_with_tax_objs[current_species].is_leaf? then 
				# nothing to do, species is already part of the tree
				next

			else
				# species not part of the tree, need the parent
				parent_node = find_parent_in_tree(current_species, lineage_root_to_species)

				if ! @taxa_with_tax_objs[parent_node].is_leaf? then 

					new_leaf = find_descendant_of_node_in_lineage(parent_node, lineage_root_to_species, all_lineages)

					if is_no_cycle_between_parent_and_node(parent_node, new_leaf) then 
						@taxa_with_tax_objs[new_leaf] = create_dummy_taxonomy_obj(new_leaf)
						link_node_to_parent(new_leaf, parent_node)
					end
					verify_links_siblings_to_parent(parent_node, new_leaf, lineage_root_to_species, all_lineages)
				end
			end
		end
	end

	def validate_tree(all_lineages)
		stack = [nil, @root]
		current_node = stack.pop
		while current_node do 
			ensure_node_is_valid(current_node, all_lineages)
			stack.push(* @taxa_with_tax_objs[current_node].descendants ) # * pushes elements of array, not array itself
			current_node = stack.pop
		end
	end

	def annotate_nodes_with_leaves_below
		current_stack = get_leaves

		next_stack = []
		while current_stack.any? do 
			current_stack.each do |current_node|

				current_node_obj = @taxa_with_tax_objs[current_node]
				if current_node_obj.is_leaf? then 
					current_node_obj.add_leaves(current_node)
				end
				current_parent = current_node_obj.ancestor
				if current_parent != "" then 
					@taxa_with_tax_objs[current_parent].add_leaves(current_node_obj.leaves_below_this_node) 
					next_stack.push(current_parent)
				end
			end
			# ancestor of root ("") is never pushed onto next_stack, thus next_stack is empty after current_stack = [root]

			current_stack = next_stack.uniq
			next_stack = []
		end
	end

	# export tree in newick format
	def export_tree(alternative_names)
# root:
# open bracket
# each desc:
# 	if no desc then 
# 		write "name:len,"
# 	end
# 	if excatly one desc then 
# 		desc: ... (same as each desc)
# 	end
# 	if more than one desc then 
# 		open bracket
# 		each desc: ...
# 		write close bracket, "name:len,"
# 	end
# write close bracket, ";"
# gsub( ",)", ")" )
		str = "("
		branch_lenght = 5

		stack = [nil, @root]
		current_node = stack.pop

		while current_node do 

			name_current_node = alternative_names[current_node] || current_node
			my_part_to_add = "#{name_current_node}:#{branch_lenght},"
			if ! @taxa_with_tax_objs[current_node].is_leaf? then 
				names_descendants = @taxa_with_tax_objs[current_node].descendants.collect do |descendant|
					alternative_names[descendant] || descendant
				end
				my_part_to_add = "(#{names_descendants.join(",")})" + my_part_to_add
			end

			if str.include?(name_current_node) then 
				str = str.gsub(name_current_node, my_part_to_add)
			else
				str += my_part_to_add
			end

			stack.push(* @taxa_with_tax_objs[current_node].descendants ) # * pushes elements of array, not array itself
			current_node = stack.pop
		end

		str += ");"
		str = str.gsub( ",,", "," )
		str = str.gsub( ",)", ")" )

		return str
	end

	def find_nodes_without_intronposition(first_taxon_with_intron, leaves_with_intron)
		nodes_without_intron = []

		stack = [nil]
		current_node = first_taxon_with_intron

		while current_node do
			leaves_below_current_node = @taxa_with_tax_objs[current_node].leaves_below_this_node
			if leaves_below_current_node.is_subset?(leaves_with_intron) then 
				# nothing to do
			elsif leaves_below_current_node.is_disjoint_set?(leaves_with_intron) then 
				nodes_without_intron.push( current_node )
			elsif leaves_below_current_node.is_overlapping_set?(leaves_with_intron) then 
				stack.push(* @taxa_with_tax_objs[current_node].descendants )
			end
			current_node = stack.pop
		end

		return nodes_without_intron
	end

	# def get_children_of_node(node)
	# 	if @taxa_with_tax_objs[node] then 
	# 		return @taxa_with_tax_objs[node].descendants
	# 	else
	# 		return []
	# 	end
	# end


	# HELPER METHODS

	def find_parent_in_tree(node, lineage_root_to_node)
		current_node = @root
		last_node = lineage_root_to_node[-1]
		my_place_not_found = true
		my_place = nil
		while my_place_not_found do 
			new_current_node = false
			@taxa_with_tax_objs[current_node].descendants.each do |current_child|
				if current_child == last_node then 
					# last node is part of the tree, and we just found its parent!
					break
				end
				if lineage_root_to_node.include?(current_child) then 
					new_current_node = true
					if current_node == current_child then 
						# attention! we are running into a never-ending loop
						Helper.abort "Cannot create phylogenetic lineage."
					end
					current_node = current_child
					break
				end
			end
			if ! new_current_node then 
				my_place_not_found = false
				my_place = current_node
			end
		end
		return my_place
	end

	def verify_links_siblings_to_parent(parent_node, current_node, lineage_root_to_current_node, all_lineages)
		@taxa_with_tax_objs[parent_node].descendants.each do |current_child|

			next if current_child == current_node

			lineage_root_to_current_child = extract_lineage_root_to_node(current_child, all_lineages)
			parent_of_current_node_and_current_child = find_last_common_ancestor_in_lineage(
				lineage_root_to_current_node, 
				lineage_root_to_current_child
			)

			if lineage_root_to_current_child.include?(current_node) then 
				# current node is ancestor of current child

				update_link_node_to_parent(current_child, current_node) # node, new_parent
			else
				# current node is no ancestor of current child

				if parent_of_current_node_and_current_child != parent_node then 
					# current_node and current_child have another last common ancestor than parent node, use this one!

					@taxa_with_tax_objs[parent_of_current_node_and_current_child] = 
						create_dummy_taxonomy_obj(parent_of_current_node_and_current_child)

					link_node_to_parent(parent_of_current_node_and_current_child, parent_node) 

					update_link_node_to_parent(current_node, parent_of_current_node_and_current_child) # node, new parent
					update_link_node_to_parent(current_child, parent_of_current_node_and_current_child)

				end # if parent_of_current_node_and_current_child != parent_node
			end # if lineage_root_to_current_node.include?(current_child)
		end # each child
	end

	def link_node_to_parent(node, parent)
		add_descendant_to_node(parent, node) # tell parent it has a new child
		set_ancestor_of_node(node, parent) # tell child it has a parent
	end

	def unlink_node_from_parent(node, parent)
		remove_descendants_from_node(parent, node) # tell parent is one child is gone
		remove_ancestor_from_node(node)	# tell child is has no parent
	end

	def update_link_node_to_parent(node, new_parent) # node, new_parent
		old_parent = @taxa_with_tax_objs[node].ancestor

		unlink_node_from_parent(node, old_parent)
		link_node_to_parent(node, new_parent) 
	end

	def add_descendant_to_node(parent, node)
		@taxa_with_tax_objs[parent].add_descendant(node)
	end
	def set_ancestor_of_node(node, parent)
		@taxa_with_tax_objs[node].set_ancestor(parent)
	end
	def remove_descendants_from_node(parent, node)
		@taxa_with_tax_objs[parent].remove_descendant(node)
	end
	def remove_ancestor_from_node(node)
		@taxa_with_tax_objs[node].set_ancestor("")
	end

	def update_root(new_root)
		@root = new_root
		remove_ancestor_from_node(@root)
	end

	def is_no_cycle_between_parent_and_node(parent, node)
		! (parent == node)
	end

	def ensure_node_is_valid(node, all_lineages)
		
		# a valid node must be root or have an ancestor
		if @taxa_with_tax_objs[node].is_root? && node != @root then 
			Helper.warn "No ancestor of node #{node} while building tree."
			
			lineage_root_to_node = extract_lineage_root_to_node(node, all_lineages)
			parent_node = find_parent_in_tree(node, lineage_root_to_node)

			link_node_to_parent(node, parent_node)
			verify_links_siblings_to_parent(parent_node, node, lineage_root_to_node, all_lineages)

		end

		# a valid node must have zero or more than two descendants
		if ! ( @taxa_with_tax_objs[node].is_leaf? || @taxa_with_tax_objs[node].descendants.size >= 2 ) then 

			new_child = find_another_descendant_of_node(node, all_lineages)

			if new_child then 
				# a regular taxon being descendant to 'node' was found

				@taxa_with_tax_objs[new_child] = create_dummy_taxonomy_obj(new_child)
				link_node_to_parent(new_child, node)
			else
				# no additional node found in lineages, add a pseudo-node instead
				Helper.warn "Adding pseudo-node to node #{node} while building tree."

				pseudo_node = "#{node}_child"
				@taxa_with_tax_objs[pseudo_node] = create_dummy_taxonomy_obj(pseudo_node)
				link_node_to_parent(pseudo_node, node)
			end
		end
	end
	# def find_path_to_root_in_tree(node)
	# 	path = [node]
	# 	current_node = node
	# 	root_not_visisted = true
	# 	while root_not_visisted do 
	# 		current_node = @taxa_with_tax_objs[current_node].ancestor
	# 		path.unshift current_node

	# 		if current_node == @root then 
	# 			root_not_visisted = false
	# 		end
	# 	end
	# 	return path
	# end

	def create_dummy_taxonomy_obj(name)
		# Taxonomy.init(name, ancestor, descendant)
		return Taxonomy.new(name, "", "")
	end

	# returns lineage from root to node (including node)
	def extract_lineage_root_to_node(node, all_lineages)
		all_lineages.each do |lineage|
			if ind_node = lineage.index(node) then 
				return lineage[0..ind_node]
			end
		end
	end

	def extract_lineages_node_to_leaves(node, all_lineages)
		found_lineages = []
		all_lineages.each do |lineage|
			if ind_node = lineage.index(node) then 
				found_lineages.push lineage[ind_node..-1]
			end
		end
		return found_lineages
	end

	# if node has other children, returns a sibling
	# otherwise, returns the direct descendant of node
	def find_another_descendant_of_node(node, all_lineages)

		current_children = @taxa_with_tax_objs[node].descendants

		lineages_node_to_leaves = extract_lineages_node_to_leaves(node, all_lineages)
		lineages_node_to_leaves.each do |lineage_node_to_leaf|
			if lineage_node_to_leaf.intersection(current_children) then 
				# this lineage is already part of the tree
				next
			end
			# lineage is not part of the tree. find a sibling
			return lineage_node_to_leaf.first
		end
		return false
	end

	def find_descendant_of_node_in_lineage(parent, lineage_root_to_species, all_lineages)

		last_common_with_nearest_sibling = parent

		current_children = @taxa_with_tax_objs[parent].descendants
		current_children.each do |current_child|
			lineage_root_to_current_child = extract_lineage_root_to_node(current_child, all_lineages)

			last_common = lineage_root_to_species.intersection( lineage_root_to_current_child ).last

			if lineage_root_to_species.index(last_common) > lineage_root_to_species.index(last_common_with_nearest_sibling) then 
				# found a closer sibling. update the last common node with the _nearest_ sibling
				last_common_with_nearest_sibling = last_common
			end
		end

		ind_first_uniq = lineage_root_to_species.index(last_common_with_nearest_sibling) + 1 # +1: first element after the last common

		return lineage_root_to_species[ind_first_uniq] 
	end

	# retruns the last common ancestor of two lineages
	def find_last_common_ancestor_in_lineage(lineage_root_to_node1, lineage_root_to_node2)
		return lineage_root_to_node1.intersection(lineage_root_to_node2).last
	end

end
