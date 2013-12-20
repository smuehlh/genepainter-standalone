# reads in a yaml file and returns a gene object
class YamlToGene

	def initialize(data, name)
		@gene = Gene.new(name) # gene object referring to all exon and intron objects
		@contigs = YAML.load( data ) # raw_data, which are different contigs
		ensure_contigs_format
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

	def to_gene
		# a gene object needs exons and introns (and the aligned sequence, which is added somewhere else)
		
		# method exons_original returns only exons, no alternative transkripts
		exons_original.each do |exon|
			exon_obj = Exon.new(exon["nucl_start"], exon["nucl_end"])
			@gene.exons << exon_obj
		end

		# method "introns" returns all (complete) introns and all uncertain introns!
		introns.each do |intron|
			pos, phase = intron_phase_and_position_in_protein(intron["nucl_start"])
			# TODO
			# position in case of sequence shift
			intron_obj = Intron.new(pos,
				intron["seq"].size, # length of intron
				phase # intron phase
				)
			@gene.introns << intron_obj
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

	# first nucleotide of gene has position 0, so "nucl_start" gives position in dna sequence and phase
	def intron_phase_and_position_in_protein(nucleotide_pos)
		# starting nucleotide modulus 3 equals the intron phase (3 nucleotides code for 1 amino acid)
		phase = nucleotide_pos.to_i % 3

		# starting nucleotide equals intron start
		pos = nucleotide_pos.to_i
		# pos = correct_intron_position_in_protein_for_phase(pos, phase) # only if pos in amino acids

		return pos, phase
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

	def tandem_gene_exons
	return self.exons.select{|exon| exon["tandem_gene"]}
	end

	def tandem_genes
		tandem_genes = []

		tandem_gene_number = 1
		loop do
			exons_of_tandem_gene = self.tandem_gene_exons.select{|tandem_gene_exon| tandem_gene_exon["tandem_gene"].to_i == tandem_gene_number}
			break if exons_of_tandem_gene.empty?
			tandem_genes << exons_of_tandem_gene
			tandem_gene_number += 1
		end
		return tandem_genes
	end

	def exon_length
		#    puts "DEBUG: exon_length: #{@name}: #{exons.select{|exon| exon["seq"].length != (exon["dna_end"]-exon["dna_start"]).abs}.collect{|exon| exon["number"]}.inspect}"
		exons.sum{|exon| exon["seq"] ? exon["seq"].length : 0}
	end

	## return merged ranges of the nucleotides that were found
	## e.g. [66..90, 165..204]
	def matched_nucl_ranges
		ranges = []
		ranges = exons.collect { |e| (e["nucl_start"].to_i..e["nucl_end"].to_i) if e["nucl_start"].to_i <= e["nucl_end"].to_i}
		ranges_merged = []
		ranges.each do |r|
			if ranges_merged[-1] && (m = ranges_merged[-1] + r) then
				ranges_merged[-1] = m
			else
				ranges_merged << r
			end
		end
		return ranges_merged
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

	def matched_nucls
		matched_nucl_ranges.collect do |r| r.to_a end.flatten.uniq
	end

	def unmatched_nucls
		(1..query_length_nuc).to_a - matched_nucls
	end

	def nucl_match_ratio
		(matched_nucls.length.to_f - (mismatches.length.to_f + seqshifts.length.to_f)) / query_length_nuc.to_f
	end

	def genomic_dna
		matchings.select{|e| e["type"] != "exon_alternative" || e["alternative_type"] == "mutual_exclusive" || e["alternative_type"] == "duplicated_exon"}.collect do |m| 
			m["seq"]
		end.join('')
	end
    
	def upstream_dna
		upstream_dna = @contigs.collect do |c| c["upstream"] end.compact[0]
		return upstream_dna
	end

	def downstream_dna
		downstream_dna = @contigs.collect do |c| c["downstream"] end.compact[0]
		return downstream_dna
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
  
	def gap_contig_lengths
	    lens = []
	    len = 0

	    #debugger
	    # puts"@contigs:#{@contigs.length}"
	    #leading gap?
	    e1 = @contigs[0]["matchings"].clone.sort_by{|m| m["dna_start"]}[0]
	    if e1["type"] == "exon" || e1["type"] == "exon_alternative" then
	      if e1["gen_start"] then
	        len = e1["gen_start"] - 0
	    else
	        len = e1["nucl_start"] - 0
	      end
	    end
	    lens << len

	    @contigs[0..-2].each_index do |i|
	      len = 0
	      #  puts "i:#{i}"
	      e1 = @contigs[i]["matchings"].clone.sort_by{|m| m["dna_end"]}[-1]
	      e2 = @contigs[i + 1]["matchings"].clone.sort_by{|m| m["dna_start"]}[0]
	      if (e1["type"] == "exon" || e1["type"] == "exon_alternative") && (e2["type"] == "exon" || e2["type"] == "exon_alternative") then
	        len = e2["nucl_start"] - e1["nucl_end"]
	      end
	      lens << len
	    end

	    len = 0
	    e1 = @contigs[-1]["matchings"].clone.sort_by{|m| m["dna_end"]}[-1]
	    if (e1["type"] == "exon" || e1["type"] == "exon_alternative") then
	      if @contigs[-1]["gen_len"] then
	        len = @contigs[-1]["gen_len"] - e1["gen_end"] if e1["gen_end"]
	    else
	        len = @contigs[-1]["prot_len"]*3 - e1["nucl_end"] if @contigs[-1]["prot_len"]
	      end
	    end
	    lens << len #if len > 0

	    return lens
	end

	# Injection to handle DNA stretches
	class Range
	  
	  # Returns true if the ranges overlap
	  def overlap?(other)
	    ol = false
	    #|------|
	    #   X------|
	    # if self.member?(other.begin) || self.member?(other.end) then ol = true end
	    #  |------|
	    #|------X
	    if other.member?(self.begin) || other.member?(self.end) then ol = true end
	    return ol
	  end
	  
	  # Compare by begin
	  def <=>(other)
	    self.begin <=> other.begin
	  end
	  
	  # Returns a new range if the ranges overlap or are adjacent.
	  # The resulting range is reverse if both input ranges are reverse.
	  # Returns false if the input ranges have a gap between them.
	  def +(other)
	    if self.grow(1).overlap?(other) then
	      # determin if concat is reverse
	      rev = (self.begin > self.end && other.begin > other.end)
	      b = [self.begin, self.end, other.begin, other.end].min
	      e = [self.begin, self.end, other.begin, other.end].max
	      if rev then b, e = e, b end
	      ret = (b..e)
	    else ret = false end
	    return ret
	  end

	  # Expands the edges of the range by 'by'
	  def grow(by)
	    if    self.begin < self.end then g = ((self.begin - by)..(self.end + by))
	    elsif self.begin > self.end then g = ((self.begin + by)..(self.end - by))
	    else g = ((self.begin - by)..(self.begin + by)) end
	    return g
	  end
	  
	end
end