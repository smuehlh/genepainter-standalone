# draws aligned genes
class GeneAlignment2svg

	def initialize(aligned_genes, is_default_color_scheme, is_nested_svg)
		@aligned_genes = aligned_genes
		@is_default_color_scheme = is_default_color_scheme
		@is_nested_svg = is_nested_svg

		@is_most_accurate_output = true # false: much better to understand, although the exon-streches are always at end of featuress

		prepare_genes_for_drawing # sets @all_intronpos_with_maxlength, @max_x_pos_exon, @scaling_factor_introns
	end

	# collect all intron positions and the lenght of the longest intron
	# because in the drawing, introns need to be kind of aligned
	# key: intron pos in alignment; value: length of longest intron in nt
	def prepare_genes_for_drawing

		# prepare maximal lenght -> needed to scale introns correctly
		calc_max_x_pos_exon

		if @is_default_color_scheme then
			fixed_intron_size = nil
		else
			fixed_intron_size = calc_intron_size_for_nondefault_color_scheme
		end

		# prepare intron data
		# find all intron positions in alignment (over all genes) and intron of maximal lenght
		@all_intronpos_with_maxlength = {} # key: intron position; value: lenght of longest intron

		@aligned_genes.each do |gene|

			# get list of all introns and their lenght in this gene
			intronpos_length = gene.get_all_introns_with_length

			# collect data
			intronpos_length.each do |pos_length|
				pos = pos_length[0]
				if @is_default_color_scheme then 
					len = enlongate_intron_to_avoid_gaps(pos_length[1])
				else
					len = fixed_intron_size
				end

				max_len_so_far = @all_intronpos_with_maxlength[pos]
				if max_len_so_far then
					# already found an intron at this position
					@all_intronpos_with_maxlength[pos] = len if len > max_len_so_far
				else
					# new intron
					@all_intronpos_with_maxlength[pos] = len
				end
			end

		end

		# prepare the scaling factor for introns
		calc_scaling_factor_introns

		# apply scaling factor
		@all_intronpos_with_maxlength.each do |k,v|
			@all_intronpos_with_maxlength[k] = scale_down_intron_length(v)
		end

	end

	# iterate though all genes to draw every feature
	# draw only exons, gaps in exons and introns, but not the gaps in introns
	# for these 'gaps in introns', create large box in background
	def create_svg

		svg = []

		# need max x-position for scaling: this is result of exons and introns, introns need to be scaled
		x_pos_max = calc_max_x_pos
		n_genes = @aligned_genes.size
		y_pos_inside_nested_element = 0.0

		svg_obj = Svg.new(n_genes, x_pos_max, @is_default_color_scheme) # svg object needed for actual drawing
		svg <<  svg_obj.print_header

		@aligned_genes.each_with_index do |gene, y_pos_nested_element|
			
			# index equals the y_position in the final drawing, it is used only in nested svg
			# all other elements are positions relative to position of nested svg


			if @is_nested_svg then 
			 	svg << svg_obj.print_nested_svg_header(y_pos_nested_element)
			 	y_pos = y_pos_inside_nested_element
			else
			 	y_pos = y_pos_nested_element
			end

			# draw text, crop text to maximal size
			svg << svg_obj.print_text( crop_gene_names_for_svg(gene.name), y_pos )

			# draw the background: artificially streches in exons to make the stay aligned
			background_col = svg_obj.colors[:exon_streched]
			svg << svg_obj.draw_horizontal_line(0, x_pos_max, y_pos, {type: "exon_streched"})

			# draw gene, without gaps
			gene.exons.each_with_index do |exon, ind|

				# exon: must be splitted, due to introns in other genes
				exon_startpos_drawing = calc_pos_drawing( exon.start_pos_in_aligned_protein, "exon" )
				exon_endpos_drawing = calc_pos_drawing( exon.end_pos_in_aligned_protein, "exon-ende" )
				exon_length = exon_endpos_drawing - exon_startpos_drawing

				if @is_most_accurate_output then 
					exon_pieces_startpos_drawing_with_length = 
						split_exons_at_foreign_intronpos( exon.start_pos_in_aligned_protein, exon.length_in_alignment )
					exon_pieces_startpos_drawing_with_length.each do |startpos_drawing, length_drawing|
						svg << svg_obj.draw_box( startpos_drawing, length_drawing, y_pos, svg_obj.colors[:exon] )
					end
				else
					svg << svg_obj.draw_box( exon_startpos_drawing, exon_length, y_pos, svg_obj.colors[:exon])
				end

				# intron
				intron = gene.introns[ind]
				next if intron.nil? # there is no intron after the last exon
				# offset: do not include this intron
				intron_startpos_drawing = calc_pos_drawing( intron.pos_last_aa_in_aligned_protein_before_intron, "intron" ) 
				if @is_default_color_scheme then 
					intron_length = scale_down_intron_length( enlongate_intron_to_avoid_gaps( intron.n_nucleotides ) )
					intron_color = svg_obj.colors[:intron]
				else
					intron_length = calc_intron_size_for_nondefault_color_scheme
					intron_color = determine_intron_color( intron.pos_last_aa_in_aligned_protein_before_intron )
				end

				svg << svg_obj.draw_box(intron_startpos_drawing, intron_length, y_pos, intron_color)

				# intron-gaps
				# is there an longer intron at same position?
				max_len_this_intron_pos = @all_intronpos_with_maxlength[ intron.pos_last_aa_in_aligned_protein_before_intron ]
				if max_len_this_intron_pos > intron_length then 
					gap_startpos_drawing = intron_startpos_drawing + intron_length
					gap_length = enlongate_intron_to_avoid_gaps( max_len_this_intron_pos - intron_length )

					svg << svg_obj.draw_box(gap_startpos_drawing, gap_length, y_pos, svg_obj.colors[:intron_gap], {draw_smaller_box: true})
				end

			end

			# draw gaps
			# draw them after then exons are drawn, as gaps may be on top of exons
			# split gaps because of introns in other genes
			gene.get_all_gaps_in_aligned_seq.each do |pos, len|
				gap_startpos_drawing = calc_pos_drawing( pos, "exon" )
				gap_endpos_drawing = calc_pos_drawing( pos+len, "exon-ende" )
				gap_length = gap_endpos_drawing - gap_startpos_drawing

				if @is_most_accurate_output then 
					gap_pieces_startpos_drawing_with_length = 
						split_exons_at_foreign_intronpos( pos, len )

					gap_pieces_startpos_drawing_with_length.each do |startpos_drawing, length_drawing|
						svg << svg_obj.draw_box( startpos_drawing, length_drawing, y_pos, svg_obj.colors[:exon_gap], {draw_smaller_box: true} )
					end
				else
					svg << svg_obj.draw_box( gap_startpos_drawing, gap_length, y_pos, svg_obj.colors[:exon_gap], {draw_smaller_box: true})
				end

			end

			if @is_nested_svg then 
				svg << svg_obj.print_nested_svg_footer
			end
		end

		# add legend 
		y_pos_nested_element = @aligned_genes.size + 1
		if @is_nested_svg then 
		 	svg << svg_obj.print_nested_svg_header( y_pos_nested_element )
		 	y_pos = y_pos_inside_nested_element
		else
			y_pos = y_pos_nested_element	
		end
		length_one_nt_exon = 1/3.0 # exon positions are in amino acid count
		length_one_nt_intron = 1 * @scaling_factor_introns # intron positions are in nt count, but scaled

		svg << svg_obj.print_legend(y_pos, length_one_nt_exon, length_one_nt_intron )

		if @is_nested_svg then 
			svg << svg_obj.print_nested_svg_footer
		end

		svg << svg_obj.print_footer

		return svg.join("\n")
	end

	# WARNING: all changes here might affect the webserver !!!
	def create_svg_merged_genestructure
		svg = []

		# need max x-position for scaling: this is result of exons and introns, introns need to be scaled
		x_pos_max = calc_max_x_pos
		n_genes = 1 # draw only merged structure
		y_pos = 0.0 # y_pos needs to start with 0; in this case regardless of nested or not!

		svg_obj = Svg.new(n_genes, x_pos_max, @is_default_color_scheme, true) # true: do not add extra space to canvas for legend!
		svg <<  svg_obj.print_header
		if @is_nested_svg then 
		 	svg << svg_obj.print_nested_svg_header(y_pos)
		end

		# actual drawing

		# text
		svg << svg_obj.print_text( crop_gene_names_for_svg("Merged"), y_pos )

		# background
		svg << svg_obj.draw_box( 0, x_pos_max, y_pos, svg_obj.colors[:exon], {draw_smaller_box: true, class_name: "background"})

		# introns
		@all_intronpos_with_maxlength.keys.sort.each_with_index do |pos, ind|
			intron_startpos_drawing = calc_pos_drawing( pos, "intron" )
			class_name = "intron_#{ind}"
			if @is_default_color_scheme then 
				intron_length = @all_intronpos_with_maxlength[pos]
				intron_color = svg_obj.colors[:intron]
			else
				intron_length = calc_intron_size_for_nondefault_color_scheme
				intron_color = determine_intron_color( pos )
			end

			svg << svg_obj.draw_box(intron_startpos_drawing, intron_length, y_pos, intron_color, {class_name: class_name})
		end
		if @is_nested_svg then 
			svg << svg_obj.print_nested_svg_footer
		end
		svg << svg_obj.print_footer

		return svg.join("\n")

	end

	# common introns (drawn underneath each other) will have the same color
	# thus, the position of the intron in the original alignment determines its color
	def determine_intron_color(pos)
		all_available_colors = Svg.get_list_of_intron_colors
		intron_ind = @all_intronpos_with_maxlength.keys.sort.index(pos)

		# map intron via modulo operation onto colors
		return all_available_colors[ intron_ind % all_available_colors.size ]
	end

	# introns are usually much, much longer than exons. Thus, the must be scaled down. 
	def scale_down_intron_length(len)
		return len * @scaling_factor_introns 
	end

	def calc_scaling_factor_introns
		if @is_default_color_scheme then 
			total_space = @max_x_pos_exon / Svg.get_exon_intron_ratio
			final_max_x_pos_intron = total_space * Svg.get_exon_intron_ratio
			return @scaling_factor_introns = final_max_x_pos_intron / calc_max_x_pos_intron
		else
			return @scaling_factor_introns = 1
		end
	end
	def calc_intron_size_for_nondefault_color_scheme
		return ( @max_x_pos_exon * (5.0/1000) ).round(2) # each intron should account for 0.5% of total exon length
	end
	def calc_max_x_pos
		# sum of exon and intron
		@max_x_pos_exon + calc_max_x_pos_intron
	end
	def calc_max_x_pos_intron
		if @is_default_color_scheme then 
			# introns must be drawing in proportion to their actual size
			return @all_intronpos_with_maxlength.values.sum
		else
			# introns must be drawn in fixed size
			return @all_intronpos_with_maxlength.keys.size * calc_intron_size_for_nondefault_color_scheme
		end
	end
	def calc_max_x_pos_exon
		# the last position to be drawn is the end of the right-most exon

		x_max_exons =  @aligned_genes.collect {|gene| gene.aligned_seq.size}.max # maximum space needed for exons equals the space needed for the alignment
		@max_x_pos_exon = x_max_exons
	end

	# offset at a certain position due to introns (length) occuring until this position
	def calc_offset_due_to_introns(pos, is_include_pos)
		introns_maxlen_upto_pos = @all_intronpos_with_maxlength.select { |intronpos, intronlen| intronpos < pos }
		offset_upto_pos = introns_maxlen_upto_pos.values.sum # sum up all values (the maximal lenght at this pos)

		if is_include_pos && @all_intronpos_with_maxlength.has_key?(pos) then 
			offset_upto_pos += @all_intronpos_with_maxlength[pos]
		end
		
		return offset_upto_pos
	end

	# add the offset due to introns 
	# handle exons and introns differently: 
	# 	exons, offset needs to be calculated including the pos itself
	# 	introns, offset should not include the pos itself
	def calc_pos_drawing(pos, type)
		if type == "exon" then 
			is_exon = true
		else
			is_exon = false
		end
		return pos + calc_offset_due_to_introns(pos, is_exon)
	end
	def get_intronpos_within_range(startpos, endpos)
		@all_intronpos_with_maxlength.keys.select do |key|
			if startpos <= key && key < endpos then 
				key
			end
		end.sort 
	end

	# split exon (and gaps) into pieces by introns in other genes
	# convert new positions and lengthes into drawing_positions!
	def split_exons_at_foreign_intronpos( startpos, len )
		startpos_len_pieces = []

		endpos = startpos + len # endpos regarding the underlying alignment, not svg

		this_start = startpos
		get_intronpos_within_range(startpos, endpos).each do |this_end|

			# add exon piece to results, convert to drawing positions
			startpos_drawing = calc_pos_drawing( this_start, "exon" ) # offset: include intron located at start pos
			endpos_drawing = calc_pos_drawing( this_end, "exon-ende" ) # offset: don't include intron located at end pos
			length_drawing = endpos_drawing - startpos_drawing
			startpos_len_pieces << [startpos_drawing, length_drawing]

			this_start = this_end
		end

		# add last exon piece
		startpos_drawing = calc_pos_drawing( this_start, "exon" )
		endpos_drawing = calc_pos_drawing( endpos, "exon-ende" )
		length_drawing = endpos_drawing - startpos_drawing
		startpos_len_pieces << [startpos_drawing, length_drawing]

		return startpos_len_pieces
	end

	def enlongate_intron_to_avoid_gaps(num)
		return num + 1
	end

	# crops gene names to fixed length (which is also used in all other output formats)
	def crop_gene_names_for_svg(name)
		return name[0...GeneAlignment.max_length_gene_name]
	end
end