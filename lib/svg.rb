class Svg
	# has all parameters and ratios for drawing
	# talks to the SvgPainter

	def initialize(n_lines, x_pos_max, is_default_color_scheme)
		@height_per_gene = calc_height_per_gene(n_lines)
		@scaling_factor = calc_scaling_factor_per_gene(x_pos_max)
		@font_size = @height_per_gene.floor

		@colors = is_default_color_scheme ? normal_colors : fancy_intron_colors
	end

	def calc_height_per_gene(n_lines)
		return Svg.parameters[:total_height] / n_lines * Svg.ratios[:gene_to_space]
	end

	def calc_scaling_factor_per_gene(x_pos_max)
		return (Svg.parameters[:total_width] - Svg.parameters[:distance_genename_to_gene]) / x_pos_max
	end

	def scale_and_shift_pos(pos)
		# scales position by scaling factor and shifts it by offset for gene name
		# this must be done for all position-data of boxes
		return (pos * @scaling_factor) + Svg.parameters[:distance_genename_to_gene]
	end

	def scale_lenght_specification(pos)
		# scales position by scaling factor, but does not shift by offset
		# shifting must not be done for length specification
		return pos * @scaling_factor
	end

	def shift_text_pos(pos)
		# shifts the baseline of text
		# effecitvly, the text is centered towards the baseline
		return pos + @font_size * (3.0/4.0)
	end
	def round_pos(val)
		# the svg will be full of rounding inprecisions anyway, but with its values rounded, its at least human readable
		return val.round(4)
	end

	# draws boxes
	# gets the 'real' positions and scales them down to the actual viewbox
	def draw_box(x_start, x_length, y_pos, type)
		x_start_draw = round_pos( scale_and_shift_pos(x_start) )
		x_len_draw = round_pos( scale_lenght_specification(x_length) )
		y_start_draw = round_pos( y_pos * @height_per_gene )
		y_len_draw = round_pos( @height_per_gene * Svg.ratios[:gene_to_space] )

		case type
		when "exon"
			this_col = @colors[:exon_default]
		when "exon-gap"
			this_col = @colors[:exon_extension]
		when "intron" 
			this_col = @colors[:intron_default]
			if this_col.kind_of?(Array) then 
# debugger
puts "welche farbe???"
				this_col = this_col.first
			end
		else 
			this_col = @colors[:intron_extension]
		end
		return Painter.box( x_start_draw, y_start_draw, x_len_draw, y_len_draw, this_col )
	end

	def print_text(txt, y_pos)
		y_pos_draw = round_pos( shift_text_pos( y_pos * @height_per_gene) )
		return Painter.text(Svg.parameters[:x_pos_genename], y_pos_draw, @font_size, txt)
	end

	def normal_colors
		# standard color scheme
		return {
			exon_default: '#e66101', # "darkorange",
			exon_extension: '#fdb863', # "sandybrown",
			intron_default: '#5e3c99', # "mediumpurple",
			intron_extension: '#b2abd2' # "thistle"
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
			intron_extension: "white"
		}
	end
	def self.parameters
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
	def self.ratios
		# definition of ratios between elements in viewbox
		return {
			gene_to_space: 0.8, # ensure some space between two "lines" of genes
			exons_to_introns: 0.5, # exons and introns not drawn to scale but such that they fit into the gene images
			legend_to_gene: 3.0 # legends as large as three gene images
		}
	end

	def self.get_exon_intron_ratio
		return ratios[:exons_to_introns]
	end
	def self.get_viewbox_size
		return { width: parameters[:total_width], height: parameters[:total_height] }
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
				">"\
				"#{txt}"\
				"</text>"
				# "dominant-baseline=\"middle\"" \
		end

		def header(width=1000,height=700, viewbox_width=Svg.get_viewbox_size[:width], viewbox_height=Svg.get_viewbox_size[:height])
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