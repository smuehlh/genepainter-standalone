# reads a GFF file and returns a gene object
class GffToGene
	def initialize(data, name, alignment_seq)
		@gene = Gene.new(name, alignment_seq) # gene object referring to all exon and intron objects
		@data = extract_cds_parts_from_gff(data) # reduce raw data to "important" parts
		
		# no need to save alignment_seq, as the queryseq is not included in standard GFF ...
		# -> not possbile to check if both match
	end

	# extracts all lines describing one (or the onliest) transkript
	# sets the start position of very first exon (this is necessary to have 0-based @gene objects)
	# WARNING: all changes here might affect the web server
	def extract_cds_parts_from_gff(gff)
		gff = gff.lines
		if gff.any? {|line| line.include?("ScipioResult")} then
			# gff is obtained from Scipio (and does not follow standard format)
			cds_lines = gff.select {|line| line.match("protein_match")}
		else
			first_mRNA_line = gff.find {|line| line.match("mRNA")}
			first_mRNA_id = nil # default: use every CDS description line
			if first_mRNA_line then 
				first_mRNA_id = get_id_from_attributes(first_mRNA_line)
			end
			cds_lines = gff.select do |line|
				line.match("CDS") && is_child_of(line, first_mRNA_id)
			end
		end

		if cds_lines.empty? then 
			Helper.abort "cannot parse gff."
		end
		@start_first_exon = parse_start_pos_from_gff_line(cds_lines.first)
		return cds_lines
	end

	def to_gene
		# data contains only CDS lines

		last_exon_end_cdna = 0 # end position of last exon in cdna seq ; cdna pos of first exon should be 0
		last_exon_end_gene = 0 # end position in gene (as reported in gff) of last exon
		@data.each do |line|
			start_pos, stop_pos, phase = parse_gff_line(line.chomp)

			exon_start_cdna = last_exon_end_cdna
			exon_end_cdna = exon_start_cdna + (stop_pos - start_pos + 1) # start + length

			exon_obj = Exon.new(exon_start_cdna, exon_end_cdna)
			@gene.exons.push(exon_obj)

			if last_exon_end_cdna != 0 then
				# not the very first exon => this exon has a preceeding intron
				intron_length = (last_exon_end_gene - start_pos).abs
				intron_obj = Intron.new(last_exon_end_cdna, # last nt before intron
					intron_length, # intron length
					phase # phase
					)
				@gene.introns.push(intron_obj)

			end
			last_exon_end_cdna = exon_end_cdna
			last_exon_end_gene = stop_pos
		end

		return @gene
	end

	# extract start/stop and phase from gff line
	def parse_gff_line(line)
		parts = line.split(/\t/)

		start_pos = gff2ruby_counting(parts[3])
		stop_pos = gff2ruby_counting(parts[4])

		# switch start and stop if gene is located on minus-strand
		# gff also contains strand, but scipio already switches start and stop pos if its located on the minus strand
		if start_pos > stop_pos then 
			tmp = start_pos
			start_pos = stop_pos
			stop_pos = tmp
		end

		phase = parts[7] || "?"
		# gff phase =! reading frame
		if phase == "1" then 
			phase = "2"
		elsif phase == "2"
			phase = "1"
		end	

		return start_pos, stop_pos, phase
	end
	def parse_start_pos_from_gff_line(line)
		line.split("\t")[3].to_i - 1
	end

	# start/stop values in GFF 3.0 (Standard) are 1-based
	# start/stop values in GFF 3.0 (Scipio format) are based on gene coordinates
	def gff2ruby_counting(num)
		num.to_i - 1# - @start_first_exon
	end

	def get_id_from_attributes(str)
		match_data = str.match(/ID=([^;]+);/)
		if match_data then 
			return match_data[1]
		else
			return nil
		end
	end

	# returns true, if parent_id is part of attributes list _or_ if parent_id is nil
	# returns false otherwise
	def is_child_of(attributes, parent_id)
		match_data = /parent=(#{parent_id})/i.match(attributes)
		return true if parent_id.nil?
		if match_data then 
			return match_data[1]
		else
			return false
		end
	end

end