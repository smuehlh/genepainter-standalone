module Svg
	# all the constants needed to draw an SVG 
	# and all the methods doing the printing itself
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