# reads in a yaml file and returns a gene object
class YamlToGene

	def initialize(data, name)
		@gene = Gene.new(name) # gene object referring to all exon and intron objects
		@contigs = YAML.load( data ) # raw_data, which are different contigs
		ensure_contigs_format
		ensure_exon_intron_numbers
	end

	# expect @contig to be an array containg the gene structure
	# but might also be an hash containing gene structures for more than one gene
	# in this case: extract the structure with has same name as @gene.name
	# if this is not included in @contig, abort
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

	def to_gene
		# a gene object needs exons and introns (and the aligned sequence, which is added somewhere else)

		exon_offset_due_to_seqshifts = 0 # offset due to seqshifts etc.

		# method exons_original returns only exons, no alternative transkripts
		exons_original.each do |exon|
			start_pos, stop_pos = exon["nucl_start"], exon["nucl_end"]

			exon_start_cdna = start_pos + exon_offset_due_to_seqshifts
			exon_end_cdna = exon_start_cdna + (stop_pos - start_pos) # start + length

			undetermined_pos = exon["undeterminedlist"]
			inframe_stop_pos = exon["inframe_stopcodons"]
	
			exon["seqshifts"].each do |seqshift|
	
				additional_target_seq = (seqshift["dna_end"]-seqshift["dna_start"]).abs 
				# when altering these conditions, please verify inframe-stopcodons and sequenceshifts are still recognized.
				# test with e.g. OtgMhc15 and TnMyo1Ea
				if ( seqshift["nucl_start"] == seqshift["nucl_end"] && 
						! is_phase_2(additional_target_seq) &&
						! undetermined_pos.include?(seqshift["prot_start"]) && 
						( seqshift["prot_start"] != seqshift["prot_end"] ||	
							inframe_stop_pos.include?(seqshift["prot_start"])
						)
					) || 
					( seqshift["nucl_start"] != seqshift["nucl_end"] && is_phase_1(additional_target_seq) ) then
					# genomic dna contains stop codon, that is not translated or a sequence shift.
					missing_target_seq = fill_up_codons(additional_target_seq)
					exon_end_cdna += missing_target_seq

					exon_offset_due_to_seqshifts += missing_target_seq
				end

				if (seqshift["nucl_start"] < seqshift["nucl_end"] &&
						seqshift["dna_start"] == seqshift["dna_end"]
				) then 
					# genomic dna contains additional bases, that are not translated.
					superfluous_target_seq = (exon["seqshifts"][0]["nucl_start"] - exon["seqshifts"][0]["nucl_end"]).abs
					exon_end_cdna -= superfluous_target_seq

					exon_offset_due_to_seqshifts -= superfluous_target_seq
				end

			end

			exon_obj = Exon.new(exon_start_cdna, exon_end_cdna)
			@gene.exons.push(exon_obj)

			intron = intron_by_intronnumber(exon["number"])
			# due to webscipio "Gap", intron might not exist
			if ! intron.nil? then  
			
				intron_length = intron["seq"].size
				pos, phase = intron_phase_and_position_in_dna(exon_end_cdna) # last nt before intron	

				intron_obj = Intron.new(exon_end_cdna, # last nt before intron
					intron_length, # intron length
					phase # phase
					)
				@gene.introns.push(intron_obj)
			end
		end

		return @gene
	end

	def correct_intron_position_in_protein_for_phase(pos, phase)
		# only necessary if the position is in amino acids
		if phase == 2 then
			pos += 1
		end
		return pos
	end

	# first nucleotide of gene has position 0
	# =>  "nucl_start" gives position in dna sequnce (and phase)
	def intron_phase_and_position_in_dna(nucleotide_pos)
		pos = nucleotide_pos.to_i
		phase = pos % 3
		return pos, phase
	end

	# first nucleotide of gene has position 0, so "nucl_start" gives position in dna sequence and phase
	def intron_phase_and_position_in_protein(nucleotide_pos)
		# starting nucleotide modulus 3 equals the intron phase (3 nucleotides code for 1 amino acid)
		phase = nucleotide_pos.to_i % 3

		# starting nucleotide equals intron start
		pos = nucleotide_pos.to_i
		pos = correct_intron_position_in_protein_for_phase(pos, phase) # only if pos in amino acids

		return pos, phase
	end
	
	def intron_by_intronnumber(num)
		# method "introns" returns all (complete) introns and all uncertain introns!
		introns.each do |intron|
			if intron["number"] == num then 
				return intron
			end
		end
		# no intron with corresponding number found
		return nil
	end
	def fill_up_codons(num)
		if is_phase_1(num) then 
			num += 2
		elsif is_phase_2(num)
			num += 1
		end
		# phase 3 : nothing to do, already a valid codon
		return num
	end
	def is_phase_1(num)
		n_codons = num / 3
		return n_codons * 3 + 1 == num
	end
	def is_phase_2(num)
		n_codons = num / 3
		return n_codons * 3 + 2 == num 
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