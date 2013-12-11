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

	def scale_genes_to_canvas_old(gene, factor_exons, factor_introns)
		scaled = gene.collect.with_index do |pos_length, ind|
			scaled_pos_length = []
			if ind.even? then
				# its a exon
				scaled_pos_length = pos_length.map { |ele| ele * factor_exons }
			else
				# its a intron
				scaled_pos_length = pos_length.map { |ele| ele * factor_introns }
			end
			scaled_pos_length
		end
		return scaled
	end

	def scale_genes_to_canvas(unscaled, max_exons, max_introns)
		# maximal width for a gene and how its distributed to exons and introns
		width_for_gene = parameters[:total_width] - parameters[:distance_genename_to_gene]
		width_exons = ratios[:exons_to_introns] * width_for_gene
		width_introns = (1 - ratios[:exons_to_introns]) * width_for_gene
		# simple cross-multiplication...
		scaling_factor_exons = width_exons / max_exons
		scaling_factor_introns = width_introns / max_introns 

		scaled = unscaled.each do |gene|
			gene[:exons] = gene[:exons].multiply_by(scaling_factor_exons)
			gene[:introns] = gene[:introns].multiply_by(scaling_factor_introns)
		end

		return scaled
	end

	def print_names_and_genes(names, genes, y_pos, is_many_intron_colors)
		svg_out = []
		names.each_with_index do |name, ind|
			# index fits to genes and y_pos
			gene = genes[ind]
			this_y_pos = y_pos[ind]

			svg_out << Painter.text(parameters[:x_pos_genename], this_y_pos, name)
			# svg_out << Painter.box
			# TODO
			debugger
			puts "convert gene[:exons] to boxes"
			puts "order does not matter, its ok to iterate first through all exons and then through all introns"
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
		output << Painter.text(x_pos_start, y_pos, "Legend")
		y_pos += size_per_gene
		output << Painter.text(x_pos_start, y_pos, "Exon")
		output << Painter.box(x_pos_start, y_pos, 10, height_of_gene_feature, "darkcyan")
		output << Painter.text(x_pos_middle, y_pos, "Gap in aligned exon")

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

		def text(x,y,txt)
			# dominant-baseline allows to align text vertically to line representing the gene
			return "<text x=\"#{x}\" "\
				"y=\"#{y}\" "\
				"font-size=\"10\" "\
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
				"viewBox=\"-1 -10 #{viewbox_width} #{viewbox_height}\" "\
				"preserveAspectRatio=\"xMinYMin meet\" "\
				"xmlns=\"http://www.w3.org/2000/svg\">"
		end
		def footer
			return "</svg>"
		end
	end

end