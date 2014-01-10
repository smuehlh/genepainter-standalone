class GeneAlignment2pdb

	attr_reader :pdb

	# get needed information from the alignment with gene structures
	# and handle pdb object

	Min_alignment_score = 0.7

	def initialize( alignment_genestruct, ref_prot, force_alignment, ref_struct_only, consensus, 
			path_to_pdb, pdb_chain, 
			nw_option )
		@alignment_genestruct = alignment_genestruct
		@ref_seq = extract_ref_seq_from_alignment(ref_prot) # reference sequence 
		@force_alignment = force_alignment || false # map introns regardless the alignment score
		@ref_struct_only = ref_struct_only || false  # map only introns occuring in the reference structure
		@min_conservation_val = consensus || 0.8 # map only introns conserved in at least X % of all genes

		# pdb object
		@pdb = instantiate_pdb_obj(path_to_pdb, pdb_chain)

		# needleman wunsch object
		@penalize_endgaps = nw_option || false
		@nw = nil
	end

	def extract_ref_seq_from_alignment(refseq_name)
		names, seqs_and_patterns = Sequence.convert_fasta_array_back_to_arrays( @alignment_genestruct )
		seq = ""

		# find sequence 

		if refseq_name then
			# search for ref_ seq
			if ! refseq_name.start_with?(">") then
				refseq_name = ">" + refseq_name
			end
			ind = names.index(refseq_name)
			Helper.abort("#{refseq_name} not found in alignment.") if ind.nil?
			seq = seqs_and_patterns[ind]
		else
			# return first seq in alignment
			# make sure it is not a gene structure, but a sequence
			ind = names.index {|name| name != "structure"}
			seq = seqs_and_patterns[ind]
		end

		# validate sequence

		if ! Sequence.is_protein_sequence_valid(seq) then
			Helper.abort("Cannot map gene structure onto protein structure: Invalid characters in reference sequence #{names[ind]}")
		end
		return seq
	end

	def align_pdb_seq_with_ref_seq
		# instance is aligned
		@nw = NeedlemanWunsch.new(@ref_seq, @pdb.seq, @penalize_endgaps)
		return @nw.score
	end

	def is_alignment_good_enough(score)
		return ( score >= Min_alignment_score || @force_alignment )
	end

	def map_genestruct_onto_pdb
		all_pos_frames = []
		if @ref_struct_only then
			# find the genestructure of the reference sequence
			# TODO

		else
			# get the consensus pattern
			stats = GeneAlignment::Statistics.new(@alignment_genestruct, true ) # true: input is fasta-formatted
			# TODO
		end

		#reduce pattern to introns conserved enough

	end

	def write_pymol_files
		# TODO
		# maybe use heredoc instead of template?
		# check if seq is already loaded and load it only if not?s

	end 

	def instantiate_pdb_obj(path, chain)
		data = File.readlines(path) 
		return Pdb.new(data, chain)
	end

	# template stuff here

end