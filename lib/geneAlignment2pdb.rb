class GeneAlignment2pdb

	# aligns a reference sequence from multiple sequence alignment with PDB sequence
	# maps gene structure onto the PDB sequence

	def initialize(aligned_seq, aligned_struct, options)

		@ref_seq = aligned_seq # aligned with input alignment
		@ref_structure = aligned_struct # algined with input alignment

		@pdb = init_pdb_obj( options[:path_to_pdb], options[:pdb_chain] )

		validate_sequences()
		@aligned_ref, @aligned_pdb = align_ref_seq_with_pdb_seq() # aligned with each other
	end

	# creates a PDB object; chain may be nil
	def init_pdb_obj(path_to_file, chain)
		data = File.readlines(path_to_file) # file existance already checked by OptParser
		Pdb.new( data, chain)
	end

	# make sure ref and pdb sequences contain no special characters
	# because any special chars are not part of scoring matrix - specification
	def validate_sequences
		err_msg = "Invalid characters in "
		if ! Sequence.is_protein_sequence_valid(@ref_seq) then
			err_msg << "PDB sequence. Replacing them by 'X'"
			Helper.log err_msg
			@ref_seq = Sequence.replace_invalid_chars_in_protein_sequence(@ref_seq)
		end
		if ! Sequence.is_protein_sequence_valid(@pdb.seq) then 
			err_msg << "protein sequence. Replacing them by 'X'"
			Helper.log err_msg
			@pdb.seq = Sequence.replace_invalid_chars_in_protein_sequence(@pdb.seq)
		end

	end

	def align_ref_seq_with_pdb_seq

		# a scoring object; values: match, mismatch [both useless in case of Blosum Matrix], gappenaltiy
		scoring_obj = Align::Blosum62.new(nil,nil,-5)
		# a aligner object; values: the seqs, the scoring object and the symobol used for an gap
		ref_seq_without_gaps = @ref_seq.delete("-")
		align_obj = Align::NeedlemanWunsch.new(ref_seq_without_gaps, @pdb.seq, {scoring: scoring_obj, skip_obj: "-"} )

		aligned_ref, aligned_pdb = align_obj.align

		return aligned_ref.join(""), aligned_pdb.join("")
	end

	def map_genestructure_onto_pdb

		mapped_intron_pos = []
		
		intron_regex = Regexp.new( "[0|1|2|?]" )
		all_intron_pos = Helper.find_each_index(@ref_structure, intron_regex)

		all_intron_pos.each do |pos_aligned_gene|

			if @ref_seq[pos_aligned_gene] == "-" then 
				# gap in reference sequence at intron position, no mapping possible
				Helper.log "Intron position #{pos_aligned_gene} in cannot be mapped onto PDB sequence."
				next
			end
			# 1) pos in aligned_seq (input alignment) to pos in unaligned seq
			pos_gene = Sequence.alignment_pos2sequence_pos(pos_aligned_gene, @ref_seq)
			# 2) pos in unaligned seq to pos in seq aligned with pdb
			pos_aligned_pdb = Sequence.sequence_pos2alignment_pos(pos_gene, @aligned_ref)
			if @aligned_pdb[pos_aligned_pdb] == "-" then 
				# gap in aligned pdb sequence, no mapping possible
				Helper.log "Intron position #{pos_aligned_gene} cannot be mapped onto PDB sequence."
				next
			end
			# 3) shift index + 1 as pdb indices start counting with one (ruby2human_counting)
			pos_pdb = pos_alignedpdb2pymol_pos_pdb ( pos_aligned_pdb )

			Helper.log "Intron at position #{Helper.ruby2human_counting(pos_aligned_gene)} => PDB position #{pos_pdb}"

			intron_phase = @ref_structure[pos_aligned_gene]

			mapped_intron_pos << "#{pos_pdb}:#{intron_phase}"

		end
		# in PDB, but no information about position in reference
		unmapped_aa_pos_in_aligned_pdb = Helper.find_each_index(@aligned_ref, "-")
		unmapped_aa_pos_in_pdb = unmapped_aa_pos_in_aligned_pdb.map do |pos|
			pos_alignedpdb2pymol_pos_pdb( pos )
		end

		# write mapped positions into python template
		output_color_exons = fill_color_exons_template(mapped_intron_pos, unmapped_aa_pos_in_pdb)
		output_color_splicesites = fill_color_splicesites_template(mapped_intron_pos, unmapped_aa_pos_in_pdb)

		return output_color_exons, output_color_splicesites

	end

	def fill_color_exons_template(mapped_intron_pos, unmapped_aminoacid_pos)

		# find the last mapped position
		# to color the exon after the last mapped intron
		pos_last_mapped_aa_aligned_pdb = @aligned_ref.rindex(/[^-]/) # last non-gap
		pos_last_mapped_aa_pdb = pos_alignedpdb2pymol_pos_pdb( pos_last_mapped_aa_aligned_pdb )

		faked_phase = '\'x\'' # phase does not matter here
		mapped_intron_pos = mapped_intron_pos.map{|pos_phase| pos_phase.sub(/:./, ":#{faked_phase}")}

		all_mapped_pos = mapped_intron_pos + ["#{pos_last_mapped_aa_pdb}:#{faked_phase}"]
		return Template.generate_color_exons_script( all_mapped_pos.join(", "), @pdb.chain, unmapped_aminoacid_pos.join(", ") )
	end

	def fill_color_splicesites_template(mapped_intron_pos, unmapped_aminoacid_pos)
		# escape phase if it is not a number
		mapped_pos_with_successor =
			mapped_intron_pos.collect do |pos_phase|
				pos, phase = pos_phase.split(":")
				successor_pos = pos.to_i + 1

				if Intron.is_phase_valid_but_non_descriptive(phase) then
					pos_phase = "#{pos}:\"#{phase}\""
					successor_phase = Template.successor_phase_for_nondescriptive_intron
				else
					successor_phase = phase.to_i + 10
				end
				[pos_phase, "#{successor_pos}:#{successor_phase}"]
			end
		return Template.generate_color_splicesites_script( mapped_pos_with_successor.join(", "), @pdb.chain, unmapped_aminoacid_pos.join(", ") )
	end

	# map pos in aligned seq onto unaligned seq
	# and shift pos by 1 as pymol starts counting with 1, but ruby with 0
	def pos_alignedpdb2pymol_pos_pdb(aligned_pos)
		return Helper.ruby2human_counting(
			Sequence.alignment_pos2sequence_pos(aligned_pos, @aligned_pdb)
			)
	end


	def log_alignment(ref_name, pdb_file)
		Helper.log "Underlying alignment"
		Helper.log ref_name
		Helper.log @aligned_ref
		Helper.log File.basename(pdb_file)
		Helper.log @aligned_pdb
	end
	def self.log_howtouse_pythonscripts
		Helper.log "Execute file 'color_exons.py' in PyMol to color exons in alternating colors."
		Helper.log "Execute file 'color_splicesites.py in PyMol to color only residues directly before and after introns."
		Helper.log "How to execute script 'script_name' in PyMol:"
		Helper.log "Open PDB file"
		Helper.log "[command line] run path_to/script_name"
		Helper.log "[command line] script_name"
		Helper.log "" 
	end
end