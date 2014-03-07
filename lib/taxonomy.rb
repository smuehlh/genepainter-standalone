class Taxonomy

	# a species has itself as child and is last common ancestor of itself
	# species has empty descendants-list
	# this is to assure that genes specific for a species can be assigned to this species like to an "real"(=inner node) last common ancestor
	attr_reader :name, 
		:ancestor, :descendants, 
		:distance_to_root, :children, :last_common_ancestor_of

	# input
	# name - taxon name
	# ancestor - ancestor name 
	# descendant - descendant name (descendants is list); empty if species
	# distance_to_root - [Int] number of nodes from here to root
	# child - [String] added to children list (=leaves below to this node) ; if taxon is an species, its its own child
	def initialize(name, ancestor, descendant, distance_to_root, child)
		@name = name
		@ancestor = ancestor
		@descendants = init_descendants(descendant)
		@distance_to_root = distance_to_root
		@children = [child]

		@last_common_ancestor_of = init_last_common_ancestor
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
		if is_species? then
			return [ @name ]
		else
			return []
		end
	end

	def add_child(child)
		@children |= [ child ]
	end

	def add_lca(child)
		@last_common_ancestor_of |= [ child ]
	end

	def add_descendant(child)
		@descendants |= [ child ] if ! child.empty?
	end

	def is_species?
		@descendants.empty?
	end

	def is_last_common_ancestor?
		@last_common_ancestor_of.any?
	end

end