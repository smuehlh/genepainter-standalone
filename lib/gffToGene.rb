# reads a GFF file and returns a gene object
class GffToGene
	def initialize(data, name)
		@gene = Gene.new(name) # gene object referring to all exon and intron objects
		@contigs = YAML.load( data ) # raw_data, which are different contigs
	end

	def to_gene
	end

	# TODO
	# get nucl_start ... for exons and introns
	# won't get dna or protein sequence
	# won't get lenght of protein
	# maybe no exons and introns but only "CDS", or "protein_match"
end
