# reads in a yaml file and returns a gene object
class YamlToGene

	def initialize(data, name, alignment_seq)
		@gene = Gene.new(name, alignment_seq) # gene object referring to all exon and intron objects
		@contigs = YAML.load( data ) # raw_data, which are different contigs
		@alignment_seq = alignment_seq

		ensure_contigs_format
		ensure_exon_intron_numbers
	end

	# expect @contig to be an array containg the gene structure
	# but might also be an hash containing gene structures for more than one gene
	# in this case: extract the structure with has same name as @gene.name
	# if this is not included in @contig, abort
	# NOTE: adapat is_yaml method of webserver when changing this method!
	def ensure_contigs_format
		if @contigs.kind_of?(Hash) then
			if @contigs[@gene.name] then
				# replace collection of gene structures by the one of interest
				@contigs = @contigs[@gene.name]
			elsif @contigs["ScipioResult"]
				# this yaml was downloaded by WebScipio, but is perfectly valid
				@contigs = @contigs["ScipioResult"]
			else
				# no gene structure for the expected gene
				Helper.abort "Missing gene structure for #{@gene.name} in respective yaml file."
			end
		end
	end

	def ensure_exon_intron_numbers
		exon_number = 0
		matchings.each do |match|
			if match["type"] == "exon" then
				if ! match.has_key?("number") then 
					exon_number += 1
					match["number"] = exon_number.to_i
				end
			elsif match["type"] == "intron" ||  match["type"] == "intron?" then
				if ! match.has_key?("number") then 
					match["number"] = exon_number.to_i
				end
			end
		end
	end

	# reject genes with queryseq other than alignment-sequence
	# NOTE: adapat is_yaml method of webserver when changing this method!
	def does_aligned_sequence_match
		query_seq = get_queryseq

		# mimic formatting of queryseq done by scipio:
		scipio_formatted_alignment_seq = @alignment_seq.upcase
		scipio_formatted_alignment_seq = scipio_formatted_alignment_seq.gsub("*", "X")
		scipio_formatted_alignment_seq = scipio_formatted_alignment_seq.gsub(/[^ACDEFGHIKLMNPQRSTVWYX]/, "")

		return scipio_formatted_alignment_seq == query_seq 
	end

	def to_gene
		# reject genes with queryseq other than alignment-sequence
		query_seq = get_queryseq
		if ! does_aligned_sequence_match then 
			# rejecting gene...
	    	Helper.log "#{@gene.name}: Aligned sequence does not match gene structure."

	    	throw(:error)
		end

		exon_offset_due_to_splitcodon = 0
		exons_original.each do |exon|

			# determine exon lenght based on the corresponding querysequence length and the phase
			# with this, sequence shifts are automatically accounted for!

			prot_start = exon["prot_start"]
			prot_end = exon["prot_end"]
			prot = query_seq[prot_start...prot_end]
			exon_length = prot.size * 3

			phase = exon["nucl_end"] % 3
			if phase == 2 then 
				# fix yaml: split-codon phase 2 belongs to preceeding exon
				phase = -1
			end	

			# add nucleotides of splitcodon at exon-end (phase); subtract split-codons added to last exon (at exon-start; offset)
			exon_length = exon_length + phase - exon_offset_due_to_splitcodon

			exon_start_pos = exon["nucl_start"]
			exon_stop_pos = exon_start_pos + exon_length

			exon_obj = Exon.new(exon_start_pos, exon_stop_pos)
			@gene.exons.push(exon_obj)

			# nucleotides added for split codon have to be subtracted from next exon
			exon_offset_due_to_splitcodon = phase	

			intron = get_intron_by_number(exon["number"])
			# due to webscipio "Gap", intron might not exist
			if ! intron.nil? then  
			
				intron_length = intron["seq"].size
				pos, phase = intron_phase_and_position_in_dna(exon_stop_pos) # last nt before intron	

				intron_obj = Intron.new(exon_stop_pos, # last nt before intron
					intron_length, # intron length
					phase # phase
					)	
				@gene.introns.push(intron_obj)
			end

		end

		return @gene
	end

	# first nucleotide of gene has position 0
	# =>  "nucl_start" gives position in dna sequnce (and phase)
	def intron_phase_and_position_in_dna(nucleotide_pos)
		pos = nucleotide_pos.to_i
		phase = pos % 3
		return pos, phase
	end

	def get_intron_by_number(num)
		# method "introns" returns all (complete) introns and all uncertain introns!
		introns.each do |intron|
			if intron["number"] == num then 
				return intron
			end
		end
		# no intron with corresponding number found
		return nil
	end

	# NOTE: adapat is_yaml method of webserver when changing this method!
	def get_queryseq
		@contigs.collect{|contig| contig["prot_seq"]}.join
	end

	### methods from webscipio

	def <=>(other)
		self.name <=> other.name
	end

	def gaps
		matchings.select do |match|match["type"] == "gap" end
	end
	       
	def matchings
		@contigs.collect do |c| c["matchings"] end.flatten
	end

	def introns
		matchings.select do |m| m["type"] == "intron" || m["type"] == "intron?" end
	end

	def uncertain_introns
		matchings.select do |m| m["type"] == "intron?" end
	end

	def exons
		matchings.select do |m| m["type"] == "exon" || m["type"] == "exon_alternative" end
	end

	def cdna
		exons.select{|e| e["type"] != "exon_alternative" || e["alternative_type"] == "mutual_exclusive" || e["alternative_type"] == "duplicated_exon"}.collect do |e| e["seq"] end.join('').upcase
	end

	def is_tandem_gene?
		@contigs.any?{|contig| contig["tandem_gene_number"]}
	end

	def has_tandem_genes?
		@contigs.any?{|contig| contig["tandem_genes"] && !contig["tandem_genes"].empty?}
	end

	def has_alternative_exons?
		matchings.any? {|match| match["type"] == "exon_alternative"}
	end

	def has_duplicated_exons?
		matchings.any? {|m| m["type"] == "exon_alternative" && m["alternative_type"] == "duplicated_exon"}
	end

	def alternative_exons
		matchings.select {|m| m["type"] == "exon_alternative"}
	end

	def has_mutually_exclusive_exons?
		matchings.any? {|match| match["type"] == "exon_alternative" && match["alternative_type"] == "mutual_exclusive"}
	end


	def mutually_exclusive_exons
		matchings.select {|m| m["type"] == "exon_alternative" && m["alternative_type"] == "mutual_exclusive"}
	end

	def duplicated_exons
		matchings.select {|m| m["type"] == "exon_alternative" && m["alternative_type"] == "duplicated_exon"}
	end

	def exons_original
		exons.select{|exon| exon["type"] == "exon"}
	end

	def exons_have_alternatives
		exons_to_return = []
		exons_original.each_with_index {|exon, i| exons_to_return << exon if exon_has_alternatives?(i+1)}
		return exons_to_return
	end

	def exons_have_mutually_exclusives
		exons_to_return = []
		exons_original.each_with_index {|exon, i| exons_to_return << exon if exon_has_mutually_exclusives?(i+1)}
		return exons_to_return
	end

	def mismatches
		exons.collect do |e| e["mismatchlist"] end.flatten.uniq.compact
	end

	def seqshifts
		exons.collect do |e|
			if (ss = e["seqshifts"]) then
				ss
				#ss.collect do |s| (s["nucl_start"]..s["nucl_end"]).to_a end
			else
				nil
			end
		end.flatten.uniq.compact
	end

	def genomic_dna
		matchings.select{|e| e["type"] != "exon_alternative" || e["alternative_type"] == "mutual_exclusive" || e["alternative_type"] == "duplicated_exon"}.collect do |m| 
			m["seq"]
		end.join('')
	end
    
	def translation
		exons.select{|ex| ex["type"] == "exon"}.collect do |m| m["translation"] end.join('')
	end

	#returns array with uniq putative exon numbers
	def count_alternative_exons
		matchings.select do |m| m["type"] == "exon_alternative" end.collect do |m| m["exon"] end.uniq
	end

	def exon_has_alternatives?(exonnumber)
		has_alternatives = false
		# puts exonnumber
		exons.each do |ex|
		  #  p ex
		  if (ex["type"] == "exon_alternative") && (ex["exon"] == exonnumber || (ex["exon_tuple"] && ex["exon_tuple"].include?(exonnumber)))
		    has_alternatives = true
		  end
		end
		return has_alternatives

	end

	def exon_has_mutually_exclusives?(exonnumber)
		return mutually_exclusive_exons.any? do |ex| ex["exon"] == exonnumber end
	end

	def exon_has_duplicated_exons?(exonnumber)
		return duplicated_exons.any? do |ex| ex["exon"] == exonnumber end
	end

	def exon_get_alternatives(exonnumber)
		exons.find_all {|ex| (ex["type"] == "exon_alternative") && (ex["exon"] == exonnumber) }
	end

	def exon_get_mutually_exclusives(exonnumber)
		exons.find_all {|ex| (ex["type"] == "exon_alternative" && ex["alternative_type"] == "mutual_exclusive") && (ex["exon"] == exonnumber) }
	end
  
end