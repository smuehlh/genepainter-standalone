module Svg
	# all the caluation and scaling and constants needed to draw an SVG
	extend self

	def colors
		# standard color scheme
		return {
			exon_default: "darkorange",
			exon_extension: "sandybrown",
			intron_default: "mediumpurple",
			intron_extension: "thistle"
		}
	end
	def fancy_intron_colors
		# focus is on common introns
		# no need for intron_extension, as all introns are displayed by 
		return {
			exon_default: "grey",
			exon_extension: "black",
			intron_default: ["crimson", "darkred", "mediumturquoise", "teal", 
				"yellowbreen", "darkgreen", "purple", "goldenrod", 
				"saddlebrown", "cadetblue"],
			intron_extension: nil
		}
	end

	def parameters
		# definition of standard parameters

		return {
			# size definition of the viewbox
			# viewbox defines a coordinate system, all values specified later on will be relative to this coordinate system
			# viewbox is then placed into the SVG-image of the user-definded size
			total_height: 70.0, 
			total_width: 100.0,
			# in viewbox, operating with fixed numbers is fine
			distance_genename_to_gene: 12.0, # gene names will need max. 12 units, gene itself can start behind 
			x_pos_genename: 0.0, # where to start the gene names
		}
	end

	def ratios
		# definition of ratios between elements in viewbox
		return {
			gene_to_space: 0.8, # ensure some space between two "lines" of genes
			exons_to_introns: 0.5, # exons and introns not drawn to scale but such that they fit into the gene images
			legend_to_gene: 3.0 # legends as large as three gene images
		}
	end

	def height_of_gene_feature(size_per_gene)
		# height per box describing the gene features, such that they don't overlap with next "line" (gene image)
		return size_per_gene * ratios[:gene_to_space]
	end

	def size_per_gene(n_genes)
		# height per gene, such that all genes _and the legend_ (!) fit into the viewbox
		return parameters[:total_height] / (ratios[:legend_to_gene] + n_genes)
	end

	def scale_genes_to_canvas(unscaled, max_exons, max_introns)

		# maximal width for a gene and how its distributed to exons and introns
		width_for_gene = parameters[:total_width] - parameters[:distance_genename_to_gene]
		width_exons = ratios[:exons_to_introns] * width_for_gene
		width_introns = (1 - ratios[:exons_to_introns]) * width_for_gene
		# simple cross-multiplication...
		scaling_factor_exons = width_exons / max_exons
		scaling_factor_introns = width_introns / max_introns 

		scaled = unscaled.collect do |gene|
# TODO collect the longest introns of each position and shift by their length

			shift_feature_starts_by_intron_lengthes(gene)

			gene.inject({}) do |new_hash, (key, val)|
				if key == :introns then
					new_hash[key] = val.multiply_by(scaling_factor_introns)
				else
					new_hash[key] = val.multiply_by(scaling_factor_exons)
				end
				new_hash
			end
		end

		return scaled
	end

	def shift_feature_starts_by_intron_lengthes(gene)
		accumulated_intron_length = 0
		# as accumulated_intron_length starts with 0, the first exon is practically 'ignored'
		0.upto(gene[:introns].size).each do |ind|
			# shift every gene feature by the intron-lenghtes of the preceeding introns
			gene.each do |key, val|
				if val[ind] then
					this_val = val[ind]
					gene[key][ind] = [ this_val[0] + accumulated_intron_length, this_val[1] ]
				end
			end

			# store intron lenght for next iteration
			if gene[:introns][ind] then
				accumulated_intron_length +=  gene[:introns][ind][1]
			end
		end
	end

	def calc_pos_of_upper_left_corner_from_pos_of_center(center, height)
		return center - height / 2
	end

	def print_names_and_genes(names, genes, y_pos, is_many_intron_colors)
		## parameters
		height_per_gene_feature = height_of_gene_feature(size_per_gene(names.size))
		this_colors = is_many_intron_colors ? colors : fancy_intron_colors

		svg_out = []
		names.each_with_index do |name, ind|
			# index fits to genes and y_pos
			gene = genes[ind]
			this_y_pos = y_pos[ind]
			this_y_pos_correced = calc_pos_of_upper_left_corner_from_pos_of_center(this_y_pos, height_per_gene_feature)
			svg_out << Painter.text(parameters[:x_pos_genename], this_y_pos, height_per_gene_feature, name)

# TODO restructure
# this should be possible in a more general way
			svg_out << gene[:exons].collect do |pos_length|
				col = this_colors[:exon_default]
				x_pos = parameters[:distance_genename_to_gene] + pos_length[0] # this is the start position, consider offset "distance_genename_to_gene"
				x_size = pos_length[1] #this is the width of the box, no need for an offset here
				Painter.box(x_pos, this_y_pos_correced, x_size, height_per_gene_feature, col)
			end
			svg_out << gene[:gaps].collect do |pos_length|
				col = this_colors[:exon_extension]
				x_pos = parameters[:distance_genename_to_gene] + pos_length[0]
				x_size = pos_length[1]
				Painter.box(x_pos, this_y_pos_correced, x_size, height_per_gene_feature, col)
			end
			svg_out << gene[:introns].collect do |pos_length|
				col = this_colors[:intron_default]
				x_pos = parameters[:distance_genename_to_gene] + pos_length[0]
				x_size = parameters[:distance_genename_to_gene] + pos_length[1]
				Painter.box(x_pos, this_y_pos_correced, x_size, height_per_gene_feature, col)
			end
		end

		return svg_out
	end

	def print_legend(size_per_gene, height_of_gene_feature)
		output = []
		y_pos = parameters[:total_height] - (size_per_gene * ratios[:legend_to_gene])
		x_pos_start = 0
		x_pos_middle = 25.0
		x_pos_end = 50.0
		output << Painter.line(x_pos_start, y_pos, parameters[:total_width], y_pos, "black")
		y_pos += size_per_gene
		output << Painter.text(x_pos_start, y_pos, height_of_gene_feature, "Legend")
		y_pos += size_per_gene
		output << Painter.text(x_pos_start, y_pos, height_of_gene_feature, "Exon")
		output << Painter.box(x_pos_start, y_pos, 10, height_of_gene_feature, "darkcyan")
		output << Painter.text(x_pos_middle, y_pos, height_of_gene_feature, "Gap in aligned exon")

		return output
	end

	module Painter
		# actual painting of the SVG
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

		def text(x,y,font,txt)
			# dominant-baseline allows to align text vertically to line representing the gene
			return "<text x=\"#{x}\" "\
				"y=\"#{y}\" "\
				"font-size=\"#{font}\" "\
				"dominant-baseline=\"middle\"" \
				">"\
				"#{txt}"\
				"</text>"
		end

		def header(width=1000,height=700, viewbox_width=Svg.parameters[:total_width], viewbox_height=Svg.parameters[:total_height])
			return "<svg version=\"1.1\" "\
				"baseProfile=\"full\" "\
				"width=\"#{width}\" "\
				"height=\"#{height}\" "\
				"viewBox=\"0 0 #{viewbox_width} #{viewbox_height}\" "\
				"preserveAspectRatio=\"xMinYMin meet\" "\
				"xmlns=\"http://www.w3.org/2000/svg\">"
		end
		def footer
			return "</svg>"
		end
	end

end