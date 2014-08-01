class Svg
	# has all parameters and ratios for drawing
	# talks to the SvgPainter
	attr_reader :colors

	@total_height = 100
	def self.total_height=(x)
		@total_height = x
	end
	def self.total_height
		@total_height
	end

	def initialize(n_lines, x_pos_max, is_default_color_scheme, no_legend=false)
		@scaling_factor = calc_scaling_factor_per_gene(x_pos_max)
		calc_canvas_height(n_lines, no_legend) # this method sets class variable total_heigth

		@colors = is_default_color_scheme ? normal_colors : fancy_intron_colors
	end

	def calc_scaling_factor_per_gene(x_pos_max)
		return (Svg.parameters[:total_width] - Svg.parameters[:distance_genename_to_gene]) / x_pos_max
	end

	# canvas is of fixed width
	# canvas height is proportional to number of lines (+ additional lines for legend)
	def calc_canvas_height(n_lines, no_legend)
		if no_legend then 
			Svg.total_height=( Svg.parameters[:height_per_gene] * n_lines )
		else
			Svg.total_height=( Svg.parameters[:height_per_gene] * (n_lines + Svg.ratios[:legend_to_gene]) )
		end
		# Svg.total_height=( Svg.parameters[:height_per_gene] * (n_lines + Svg.ratios[:legend_to_gene]) )
	end

	def scale_and_shift_x_pos(pos)
		# scales position by scaling factor and shifts it by offset for gene name
		# this must be done for all position-data of boxes
		return round_pos( (pos * @scaling_factor) + Svg.parameters[:distance_genename_to_gene] )
	end

	def scale_x_lenght(pos)
		# scales position by scaling factor, but does not shift by offset
		# shifting must not be done for length specification
		return round_pos( pos * @scaling_factor )
	end
	def scale_y_pos(pos)
		round_pos( pos * Svg.parameters[:height_per_gene] )
	end
	def get_y_length
		round_pos( Svg.parameters[:height_per_gene] * Svg.ratios[:gene_to_space] )
	end
	def shift_text_pos(pos)
		# shifts the baseline of text
		# effecitvly, the text is centered towards the baseline
		return round_pos(pos + Svg.parameters[:font_size] * (3.0/4.0))
	end
	def round_pos(val)
		# the svg will be full of rounding inprecisions anyway, but with its values rounded, its at least human readable
		return val.round(4)
	end
	def shift_y_pos_line(pos)
		# shift y_pos (is upper left corner) towards center of a box
		round_pos( (pos * Svg.parameters[:height_per_gene]) + (Svg.parameters[:height_per_gene] - get_line_width) / 2 )
	end
	def shift_y_pos_small_box(pos, ratio)
		round_pos( pos * Svg.parameters[:height_per_gene] + ((Svg.parameters[:height_per_gene] - Svg.parameters[:height_per_gene] * ratio)/2) )
	end
	def get_line_width 
		round_pos( Svg.parameters[:height_per_gene] * 0.1 )
	end

	def print_header
		return Painter.header( Svg.parameters[:total_width], Svg.parameters[:total_height] )
	end

	def print_footer
		return Painter.footer
	end
	def print_nested_svg_header(y_pos)
		y_pos_draw = scale_y_pos(y_pos)
		return Painter.nested_svg_header(y_pos_draw)
	end
	def print_nested_svg_footer
		return Painter.footer
	end

	# draws boxes
	# gets the 'real' positions and scales them down to the actual viewbox
	# opts [Hash] {:dont_scale_x [Boolean]} specifies if x position should be scaled
	# opts [Hash] {:draw_smaller_box [Boolean]} specifies if box height should be smaller than usual 
	# opts [Hash] {:class_name [String]} set a class for this box
	def draw_box(x_start, x_length, y_pos, color, opts={})
		if opts[:dont_scale_x] then 
			x_start_draw = x_start
			x_length_draw = x_length
		else
			x_start_draw = scale_and_shift_x_pos(x_start)
			x_length_draw = scale_x_lenght( x_length )
		end
		if opts[:draw_smaller_box] then 
			# box height is smaller than usual box height
			y_start_draw = shift_y_pos_small_box(y_pos, Svg.ratios[:smaller_box_to_normal_box])
			y_length = get_y_length * Svg.ratios[:smaller_box_to_normal_box]
		else
			y_start_draw = scale_y_pos(y_pos)
			y_length = get_y_length
		end
		if opts[:class_name] then 
			class_name = opts[:class_name]
		else
			class_name = ""
		end

		return Painter.box( x_start_draw, y_start_draw, x_length_draw, y_length, color, class_name )
	end

	def print_text(txt, x_pos=Svg.parameters[:x_pos_genename], y_pos)
		y_pos_draw = shift_text_pos( y_pos * Svg.parameters[:height_per_gene])
		return Painter.text(x_pos, y_pos_draw, Svg.parameters[:font_size], txt)
	end

	# does the actual scaling of x and y positions
	def draw_line(x_start, x_stop, y_start, y_stop, line_width, color, is_scale_xpos)
		if is_scale_xpos then 
			x_start_draw = scale_and_shift_x_pos( x_start )
			x_stop_draw = scale_and_shift_x_pos( x_stop )
		else
			x_start_draw = x_start
			x_stop_draw = x_stop
		end
		y_start_draw = shift_y_pos_line( y_start )
		y_stop_draw = shift_y_pos_line(y_stop )

		return Painter.line(x_start_draw, y_start_draw, x_stop_draw, y_stop_draw, line_width, color)
	end
	# opts [Hash] {:dont_scale_x [Boolean]} specifies if x position should be scaled
	# opts [Hash] {:type [String]} key of @colors hash 
	def draw_vertical_line(x_pos, y_start, y_stop, opts={})
		is_scale_xpos = ! opts[:dont_scale_x]
		color = @colors[opts[:type]] || @colors[:exon_streched]
		line_width = get_line_width

		draw_line(x_pos, x_pos, y_start, y_stop, line_width, color, is_scale_xpos)
	end
	def draw_horizontal_line(x_start, x_stop, y_pos, opts={})
		is_scale_xpos = ! opts[:dont_scale_x]
		color = @colors[opts[:type]] || @colors[:exon_streched]
		line_width = get_line_width 

		draw_line(x_start, x_stop, y_pos, y_pos, line_width, color, is_scale_xpos)
	end

	# a scalebar looks like this: 
	# |----| X bps Exon
	# |----| Y bps Intron
	def add_scalebar_to_legend(x_pos_scalebar, y_pos, length_one_nt, text)
		scalebar = []

		# number of nucleotides equaling a 'scalebar_width' units long line
		# cross multiplication to get number of nt per 10 viewbox units and scale them
		n_nt = ((1 / length_one_nt * Svg.parameters[:scalebar_width])* 1/@scaling_factor).round 
		n_nt = Helper.convert_number_to_human_readable_string(n_nt)

		x_pos_text = x_pos_scalebar + Svg.parameters[:scalebar_width] + 0.5 # add 0.5 to add some space between scalebar and the text
		scalebar << print_text("#{n_nt} #{text}", x_pos_text, y_pos)
		scalebar << draw_horizontal_line(x_pos_scalebar, Svg.parameters[:scalebar_width] + x_pos_scalebar, y_pos, {dont_scale_x: true}) # "darkgrey")
	
		y_pos_start_vertical_lines = y_pos - Svg.parameters[:height_per_gene] * 0.02 # 0.01: the line should be centered at y_pos_draw and have a total height of 2/100 of @height_per_gene
		y_pos_stop_vertical_lines = y_pos + Svg.parameters[:height_per_gene] * 0.02

		scalebar << draw_vertical_line(
			x_pos_scalebar, y_pos_start_vertical_lines, y_pos_stop_vertical_lines, {dont_scale_x: true}
		)  # color = "darkgrey")
		scalebar << draw_vertical_line(
			x_pos_scalebar + Svg.parameters[:scalebar_width], y_pos_start_vertical_lines, y_pos_stop_vertical_lines, {dont_scale_x: true}
		)# color = "darkgrey")

		return scalebar
	end

	# a legend consists of coloured boxes and their annotation and the scalebars
	# a legends needs 4 additional lines
	def print_legend(y_pos, length_one_nt_exon, length_one_nt_intron)

		box_size = 10 # 10 viewbox units look ok
		distance_text_to_box = box_size + 0.5 # text should start after the end of the box; add 0.5 to get some additional space
		distances_between_text_columns = Svg.parameters[:distance_genename_to_gene] / 2
		x_pos_2nd_col = Svg.parameters[:distance_genename_to_gene]
		x_pos_3rd_col = x_pos_2nd_col + distances_between_text_columns + distance_text_to_box
		x_pos_4rd_col = x_pos_3rd_col + distances_between_text_columns + distance_text_to_box

		is_default_color_scheme = @colors[:intron].kind_of?(Array) ? false : true

		legend = []

		legend << draw_horizontal_line(0, Svg.parameters[:total_width], y_pos, {dont_scale_x: true}) # "black"
		y_pos += 1 # continue in next line, otherwise horizontal line overlaps with text
		legend.unshift( print_text("Legend", y_pos) ) # its important that the text "Legend" appears first in SVG

		legend << print_text("Exon", x_pos_2nd_col + distance_text_to_box, y_pos)
		legend << draw_box(x_pos_2nd_col, box_size, y_pos, @colors[:exon], {dont_scale_x: true}) 
		legend << print_text("Intron", x_pos_3rd_col + distance_text_to_box, y_pos)
		# intron colors, draw them all :-)
		if is_default_color_scheme then 
			legend << draw_box(x_pos_3rd_col, box_size, y_pos, @colors[:intron], {dont_scale_x: true})
		else
			n_colors = @colors[:intron].size
			width_per_col_and_gap = box_size / n_colors.to_f
			width_per_col = round_pos(width_per_col_and_gap * 0.8)
			@colors[:intron].each_with_index do |color, ind|
				this_x_pos = x_pos_3rd_col + ind * width_per_col_and_gap
				legend << draw_box(this_x_pos, width_per_col, y_pos, color, {dont_scale_x: true})
			end
		end

		legend << add_scalebar_to_legend(x_pos_4rd_col, y_pos, length_one_nt_exon, "bps Exon")

		y_pos += 1
		legend << print_text("Exon gap", x_pos_2nd_col + distance_text_to_box, y_pos)
		legend << draw_box(x_pos_2nd_col, box_size, y_pos, @colors[:exon_gap], {dont_scale_x: true})

		if is_default_color_scheme then 
			legend << print_text("Intron gap", x_pos_3rd_col + distance_text_to_box, y_pos)
			legend << draw_box(x_pos_3rd_col, box_size, y_pos, @colors[:intron_gap], {dont_scale_x: true})

			legend << add_scalebar_to_legend(x_pos_4rd_col, y_pos, length_one_nt_intron, "bps Intron")
		end

		return legend
	end

	def normal_colors
		# standard color scheme
		return {
			exon: '#e66101', # "darkorange",
			exon_gap: '#fdb863', # "sandybrown",
			intron: '#5e3c99', # "mediumpurple",
			intron_gap: '#b2abd2', # "thistle"
			exon_streched: 'grey'
		}
	end
	def fancy_intron_colors
		# focus is on common introns
		# no need for intron_extension, as all introns are displayed by 
		return {
			exon: "darkgrey",
			exon_gap: "lightgrey",
			intron: Svg.get_list_of_intron_colors,
			intron_gap: "white",
			exon_streched: "black"
		}
	end

	def self.parameters
		# definition of standard parameters
		return {
			# size definition of the canvas and viewbox
			total_height: Svg.total_height, 
			total_width: 1000.0,
			# in this set, where size is fixed, operating with fixed numbers is fine
			distance_genename_to_gene: 125.0, # needed for 20 chars in monospace font 
			x_pos_genename: 0.0, # where to start the gene names
			height_per_gene: 10,
			font_size: 10,
			scalebar_width: 62.5
		}
	end
	def self.ratios
		# definition of ratios between elements in viewbox
		return {
			gene_to_space: 0.8, # ensure some space between two "lines" of genes
			exons_to_introns: 0.5, # exons and introns not drawn to scale but such that they fit into the gene images
			legend_to_gene: 4.0, # legends as large as three gene images
			smaller_box_to_normal_box: 0.7 # some boxes should be smaller than others
		}
	end

	def self.get_list_of_intron_colors
		return ["tomato", "firebrick", "mediumturquoise", "teal", 
				"yellowgreen", "darkgreen", "purple", "goldenrod", 
				"saddlebrown", "cadetblue"]
	end
	def self.get_exon_intron_ratio
		return ratios[:exons_to_introns]
	end

	module Painter
		# actual painting of the SVG
		extend self

		def line(x1, y1, x2, y2, width, color)
			return "<line x1=\"#{x1}\" "\
				"y1=\"#{y1}\" "\
				"x2=\"#{x2}\" "\
				"y2=\"#{y2}\" "\
				"stroke=\"#{color}\" stroke-width=\"#{width}\" "\
				"/>"
		end

		def box(x, y, width, height, color, my_class)
			return "<rect x=\"#{x}\" "\
				"y=\"#{y}\" "\
				"width=\"#{width}\" "\
				"height=\"#{height}\" "\
				"fill=\"#{color}\" "\
				"class=\"#{my_class}\" "\
				"/>"
		end

		def text(x,y,font,txt)
			# dominant-baseline allows to align text vertically to line representing the gene
			return "<text x=\"#{x}\" "\
				"y=\"#{y}\" "\
				"font-size=\"#{font}\" "\
				"font-family=\"monospace\""\
				">"\
				"#{txt}"\
				"</text>"
				# "dominant-baseline=\"middle\"" \
		end

		def header(width, height)
			return "<svg version=\"1.1\" "\
				"baseProfile=\"full\" "\
				"width=\"#{width}\" "\
				"height=\"#{height}\" "\
				"viewBox=\"0 0 #{width} #{height}\" "\
				"preserveAspectRatio=\"xMinYMin meet\" "\
				"xmlns=\"http://www.w3.org/2000/svg\">"
		end
		def nested_svg_header(y)
			return "<svg y=\"#{y}\">"
		end
		def footer
			return "</svg>"
		end
	end

end