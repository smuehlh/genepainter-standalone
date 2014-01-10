class Pdb

    attr_reader :seq

	def initialize(pdb_data, chain)
		@chain = (chain || "a").upcase
		@seq = get_seq_from_pdb_file(pdb_data)
	end

	def get_seq_from_pdb_file(data)
		seq = []

		# reduce info to lines describing atoms
		atom_lines = data.grep(/^ATOM/)

		seq_pos_last_visited = 0 # one sequence position has many atoms, defaults to 0 as sequence position count starts with 1

		atom_lines.each do |line|

			parts = line.split(/\s+/)
            # make sure that accessing fields 3-5 of array parts will work
            next if parts.size <= 6

			# reduce atom description to correct chain
			this_chain = parts[4].upcase

			if this_chain == @chain then 
				this_seq_pos = parts[5].to_i

				# build sequence, but add each amino acid only once for every atom
				if this_seq_pos != seq_pos_last_visited then 
					if this_seq_pos > seq_pos_last_visited + 1 then 
						# add gap(s)
						seq << "-" * (this_seq_pos - seq_pos_last_visited - 1)
					end
					# add amino acid in one letter code

					seq << amino_acid_3_to_1(parts[3])
                    seq_pos_last_visited = this_seq_pos
				end
			else
				# wrong chain
                # do nothing
			end
		end

        # validate sequence

        pdb_seq = seq.join("")

        if ! Sequence.is_protein_sequence_valid(pdb_seq) then
            Helper.abort("Cannot map gene structure onto protein structure: Invalid characters in PDB sequence")
        end
        return pdb_seq
	end

    def amino_acid_3_to_1(abbr)
        # TODO alle upcase
        table = {"ALA" => "A",
            "ARG" => "R",
            "ASN" => "N",
            "ASP" => "D",
            "CYS" => "C", 
            "GLU" => "E",
            "GLN" => "Q",
            "GLY" => "G",
            "HIS" => "H",
            "ILE" => "I",
            "LEU" => "L",
            "LYS" => "K",
            "MET" => "M",
            "PHE" => "F",
            "PRO" => "P",
            "SER" => "S",
            "THR" => "T",
            "TRP" => "W",
            "TYR" => "Y",
            "VAL" => "V"}

        return table[abbr.upcase]
    rescue ArgumentError
        # non standard amino acid
        return 'X'
    end
end
