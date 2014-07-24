class Taxonomy

	attr_reader :name, 
		:ancestor, :descendants, :last_common_ancestor_of, :leaves_below_this_node

	# input
	# name - taxon name
	# ancestor - ancestor name 
	# descendant - descendant name (descendants is list); empty if species
	def initialize(name, ancestor, descendant)
		@name = name
		@ancestor = ancestor
		@descendants = init_descendants(descendant)
		@last_common_ancestor_of = init_last_common_ancestor
		@leaves_below_this_node = []
	end

	# a species should have no descendants
	def init_descendants(child)
		if child.empty? then 
			return []
		else
			return [ child ]
		end
	end

	# a species should be its own last common ancestor
	def init_last_common_ancestor
		if is_leaf? then
			return [ @name ]
		else
			return []
		end
	end

	def add_lca(child)
		@last_common_ancestor_of |= [ child ]
	end

	def add_descendant(child)
		@descendants |= [ child ] if ! child.empty?
	end

	def remove_descendant(child)
		@descendants.delete(child)
	end

	def set_ancestor(ancestor)
		@ancestor = ancestor
	end

	def add_leaves(child_list)
		@leaves_below_this_node |= Array(child_list) # Array for smart conversion to array
	end

	def is_leaf?
		@descendants.empty?
	end
	def is_last_common_ancestor_of?(species)
		@last_common_ancestor_of.include?(species)
	end
	def is_root?
		@ancestor.empty?
	end

end