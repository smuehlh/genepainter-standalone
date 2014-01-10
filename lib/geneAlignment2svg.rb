# draws aligned genes
class GeneAlignment2svg

	def initialize(aligned_genes, svg_params)
		@aligned_genes = aligned_genes
		@svg_size, @is_default_color_scheme = parse_svg_parameters(svg_params)

		prepare_genes_for_drawing # sets @all_intronpos_with_maxlength, @max_x_pos_exon, @scaling_factor_introns
	end

	# params are the parsed command line arguments
	# extract requested image size and type (normal or reduced) 
	def parse_svg_parameters(params)

		# size
		if ! params[:size] then 
			param_size = { width: nil, height: nil}
		else
			param_size = { width: params[:size][:width], height: params[:size][:height] }
		end

		# type
		if ! params[:reduced] then 
			# use the default color scheme if the parameter for the color is not set at all or if it is false
			is_def_color = true
		else 
			is_def_color = false
		end
		
		return param_size, is_def_color
	end

	# collect all intron positions and the lenght of the longest intron
	# because in the drawing, introns need to be kind of aligned
	# key: intron pos in alignment; value: length of longest intron in nt
	def prepare_genes_for_drawing

		# prepare maximal lenght -> needed to scale introns correctly
		calc_max_x_pos_exon

		# prepare intron data
		# find all intron positions in alignment (over all genes) and intron of maximal lenght
		@all_intronpos_with_maxlength = {} # key: intron position; value: lenght of longest intron

		@aligned_genes.each do |gene|

			# get list of all introns and their lenght in this gene
			intronpos_length = gene.get_all_introns_with_length(true) # true: get positions in nucleotides

			# collect data
			intronpos_length.each do |pos_length|
				pos = pos_length[0]
				len = pos_length[1]

				max_len_so_far = @all_intronpos_with_maxlength[pos]
				if max_len_so_far then
					# already found an intron at this position
					@all_intronpos_with_maxlength[pos] = len if len > max_len_so_far
				else
					# new intron
					@all_intronpos_with_maxlength[pos] = len
				end
			end
break
# TODO break nur solange SVG nur fuer 1. Gen 
		end

		# prepare the scaling factor for introns
		calc_scaling_factor_introns

		# apply scaling factor
		@all_intronpos_with_maxlength.each do |k,v|
			@all_intronpos_with_maxlength[k] = fix_intron_length(v)
		end

	end

	# iterate though all genes to draw every feature
	# draw only exons, gaps in exons and introns, but not the gaps in introns
	# for these 'gaps in introns', create large box in background
	def create_svg
		svg = []
		svg << Svg::Painter.header( @svg_size[:width], @svg_size[:height] )

		# need max x-position for scaling: this is result of exons and introns, introns need to be scaled
		x_pos_max = calc_max_x_pos
		n_genes = @aligned_genes.size

		svg_obj = Svg.new(n_genes, x_pos_max, @is_default_color_scheme) # svg object needed for actual drawing
		@aligned_genes.each_with_index do |gene, y_pos|
			# index equals the y_position in the final drawing

			# draw text
			svg << svg_obj.print_text(gene.name, y_pos)

			# draw 'background' box, only visible if there is no feature (= an "gap" in an intron ) 
			svg << svg_obj.draw_box(0, x_pos_max, y_pos, "intron-gap")

			# gaps within this gene
			@all_gappos_with_length_this_gene = gene.get_gaps_with_length_mapped_onto_dna_seq

			# draw gene
			gene.exons.each_with_index do |exon, ind|
				is_calc_offset_including_pos = true

				# exon
				exon_startpos = exon.start_pos_in_dna_seq 
				exon_startpos_drawing = calc_startpos_drawing(exon_startpos, "exon")

				# strech box to include some space for gaps in exon
				exon_length = exon.length_in_nt + calc_gap_length_upto_pos(exon_startpos, exon.end_pos_in_dna_seq)
puts "Exon: #{exon.start_pos_in_dna_seq} => #{exon_startpos_drawing}: #{exon_length}"	
# TODO 
# problem 1) exongap-boxen werden ncoh nicht richtig gezeichnet (siehe ang-coro2a)
# problem 2) position intron/laenge exon: da wo hintergrund durchscheint
				svg << svg_obj.draw_box(exon_startpos_drawing, exon_length, y_pos, "exon")

				# intron
				this_intron = gene.introns[ind]
				next if this_intron.nil? # last index of exon is not matched by an intron
				intron_startpos = this_intron.pos_last_nt_in_dna_seq_before_intron
				intron_startpos_drawing = calc_startpos_drawing(intron_startpos, "intron")

				intron_length = fix_intron_length( this_intron.n_nucleotides )
puts "Intron: #{this_intron.pos_last_nt_in_dna_seq_before_intron} => #{intron_startpos_drawing}: #{intron_length}"		
puts ""		
				svg << svg_obj.draw_box(intron_startpos_drawing, intron_length, y_pos, "intron")
			end

			# gaps
			@all_gappos_with_length_this_gene.each do |pos, len|
				gap_startpos_drawing = calc_startpos_drawing(pos, "exon")
				svg << svg_obj.draw_box(gap_startpos_drawing, len, y_pos, "exon-gap")
			end

break
		end

		svg << Svg::Painter.footer

		return svg
	end

	def fix_intron_length(len)
		# introns are usually much, much longer than exons. Thus, the must be scaled down. 
		# - or -
		# their length is not important at all, just their position
		if @is_default_color_scheme then 
			# effective size is proportional to real size
			return len * @scaling_factor_introns
		else
			# effective size is fixed; not related to real size
			return ( @max_x_pos_exon * (5.0/100) ).round(2) # each intron should account for 5% of total exon length
		end
	end

	def calc_scaling_factor_introns
		total_space = @max_x_pos_exon / Svg.get_exon_intron_ratio
		final_max_x_pos_intron = total_space * Svg.get_exon_intron_ratio

		@scaling_factor_introns = final_max_x_pos_intron / calc_max_x_pos_intron
	end

	def calc_max_x_pos
		# sum of exon and intron
		@max_x_pos_exon + calc_max_x_pos_intron
	end
	def calc_max_x_pos_intron
		@all_intronpos_with_maxlength.values.sum # already in nt count
	end
	def calc_max_x_pos_exon
		# find the last postion to be drawn is the _end_ of the last exon
		# the right-most, last exon of all genes ends at nucleotide position alignment_size * 3

		x_max_exons =  @aligned_genes.collect {|gene| gene.aligned_seq.size}.max # maximum space needed for exons equals the space needed for the alignment
		@max_x_pos_exon = x_max_exons * 3 # convert aa to nt count
	end

	def calc_offset_due_to_introns(alignmentpos, is_include_alignmentpos)
		offset = 0
		@all_intronpos_with_maxlength.each do |pos, len|
			if pos < alignmentpos then 
				offset += len
			elsif pos == alignmentpos && is_include_alignmentpos 
				offset += len
			end
		end
		return offset
	end

	def calc_offset_due_to_gaps_in_exons(startpos=0, endpos)
		# problem: muss offset nur die gaps in exons beruecksichtigen, oder auch die zwischen exons?  
		offset = 0
		@all_gappos_with_length_this_gene.each do |pos, len|
			if startpos <= pos then 
				if pos < endpos then
					offset += len
				end
			end 
		end
		return offset
	end

	def calc_startpos_drawing(startpos, type)
		if type == "exon" then 
			is_exon = true
		else
			is_exon = false
		end
		return startpos + 
			calc_offset_due_to_introns(startpos, is_exon) + 
			calc_offset_due_to_gaps_in_exons(startpos)
	end

	def calc_gap_length_upto_pos(startpos, endpos)
		matching_pos_len = @all_gappos_with_length_this_gene.select do |thispos, len| 
			len if startpos < thispos && thispos < endpos
		end
		return matching_pos_len.values.sum
	end
end