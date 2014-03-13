# !/usr/bin/env ruby

if ! RUBY_VERSION.include?("2.0") then
	Helper.abort( "You are using Ruby #{RUBY_VERSION}. For GenePainter, Ruby >= 2.0 is needed.")
end

# require .rb files in library (including all subfolders)
# TODO remove ruby-debug
require 'ruby-debug'
Dir[File.join(File.dirname(__FILE__), 'lib', '**', '*.rb')].each do |file|
	require File.absolute_path(file)
end

### parse command line arguments
args = ARGV.map {|ele| ele.dup}
OptParser.backwards_compability_check(args)
options = OptParser.parse(args)

### function defintion
# write output verbosely to file
def write_verbosely_output_to_file(f_name, output)
	puts "\t writing output to #{f_name} ... "
	# add newline to string, in case it does not end with one
	if ! output.end_with?("\n") then 
		output << "\n"
	end
	IO.write( f_name, output, :mode => "w" )
end

### open log file, and specify its automatic closure at_exit
Helper.open_log( options[:path_to_log] )
at_exit { Helper.close_log }
Helper.log "Program call: #{ARGV.join(" ")}"

### read in data, use only intersect of alignment and gene structures
if options[:tax_options].any? then 
	## taxonomy, only if neccessary
	print "Read in taxonomy ..."
	taxonomy_obj = TaxonomyToGene.new( options[:tax_options][:path_to_tax], options[:tax_options][:path_to_tax_mapping] )
	puts " done."
end

print "Read in alignment and genes ..."
## alignment
aligned_seqs_names, aligned_seqs = Sequence.read_in_alignment(options[:path_to_alignment])

## genes
gene_names = [] # a list of all gene names
Dir.glob( File.join(options[:path_to_genestruct], "*.*") ) do |file|
	f_name = File.basename(file, ".*")
	gene_names << f_name
end

## create gene_objects, they will be in same order as they are in the alignment
# initialize 
gene_objects = [] # a list of all gene objects
common_names = aligned_seqs_names & gene_names

if common_names.empty? then 
	Helper.abort "No matches between alignment and gene structures."
end

# make sure all _needed_ aligned sequences are of same length
common_aligned_seqs = Sequence.ensure_seqs_have_same_length(aligned_seqs, aligned_seqs_names, common_names)
# remove common gaps if neccessary
if options[:ignore_common_gaps] then 
	common_aligned_seqs = Sequence.remove_common_gaps(common_aligned_seqs, 
		{keep_one_common_gap_of_each_set: false, insert_common_gap_between_non_gaps: false} # remove all common gaps, and never add common gaps
	) 
end

# special case: gap before/after intron
# correct intron position (default: before gap) if there is any other sequence with and intron at same pos as gap-end
introns_before_gap_pos_gene = {}

# merge gene structure and sequence into a gene object
common_names.each_with_index do |gene_name, ind|

	# find gene structure file with belongs to this gene
	matching_files = Dir.glob( File.join(options[:path_to_genestruct], "#{gene_name}.*") ) 

	# double-check that the gene was matched un-ambiguously to gene structure files
	if matching_files.size > 1 then
		Helper.abort("Cannot parse gene structure. #{gene_name} is ambiguous file name.")
	else
		file = matching_files[0]
		f_extension = File.extname(file)
	end

	# read in file and parse gene structure
	data_obj = nil
	data = IO.read(file)
	if f_extension.casecmp(".yaml") == 0 then
		# yaml format
		data_obj = YamlToGene.new(data, gene_name)

	elsif f_extension.casecmp(".gff") == 0 then
		# gff format
		data_obj = GffToGene.new(data, gene_name)

	else
		Helper.abort("Cannot parse gene structure. Unknown file type #{file}.")
	end

	# create a gene object with this structure and sequence
	gene_obj = data_obj.to_gene # method to_gene returns a gene object containing the structure
	gene_obj.add_aligned_seq(common_aligned_seqs[ind]) # ... and aligned sequence
	gene_obj.reduce_gene_to_range(options[:range]) if options[:range].any?

	gene_objects << gene_obj

	# collect all positions of introns before gap-positions in aligned sequence
	if options[:best_intron_pos] then 
		gene_obj.get_all_gap_boundaries_preceeded_by_intron.each do |intron_pos, gap_end|
			(introns_before_gap_pos_gene[intron_pos] ||= []).push( [gene_name, gap_end] )
		end
	end

end

puts " done."

if options[:best_intron_pos] && introns_before_gap_pos_gene.any? then 
	print "Corret intron positions flanking alignment gaps ..."

	Sequence.correct_introns_flanking_gaps(introns_before_gap_pos_gene, gene_objects)
	puts " done."
end

# inform which data are not used
Helper.print_intersect_and_diff_between_alignment_and_gene(aligned_seqs_names, gene_names)
if taxonomy_obj then 

	options[:tax_options][:selected_taxa] ||= [] # set empty selected_taxa explicitly to empty array
	genes_belonging_to_selected_taxa = options[:tax_options][:selected_taxa].collect do |taxon|
		taxonomy_obj.get_genes_encoded_by_taxon(taxon)
	end.flatten.uniq

	Helper.print_intersect_and_diff_between_taxonomy_genes_and_selected_taxa( 
		common_names, taxonomy_obj.known_genes, options[:tax_options][:selected_taxa], genes_belonging_to_selected_taxa
	)
end

### align genes
puts ""
print "Aligning genes ..."

# initiate an gene alignment object 
# this calculates kind of an "master format", the exon-intron-patterns for each gene
gene_alignment_obj = GeneAlignment.new(gene_objects, 
	options[:consensus], 
	options[:merge], 
	{ genes_belonging_to_selected_taxa: genes_belonging_to_selected_taxa || [], 
		is_intron_exclusive_for_selected_taxa: options[:tax_options][:is_exclusive] || false 
	},
	options[:is_long_text_based_output]
	)
puts " done."

### prepare output for every requested format
# checking for each possible format, if it was specified, is faster than looping through all specified ones
puts "Prepare output ... "

if options[:output_format_list].include?("alignment_with_intron_phases") then 
	# this is in most cases the master format
	output = gene_alignment_obj.export_as_alignment_with_introns
	f_out = options[:path_to_output] + ".fas"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("txt_simple") then
	output = gene_alignment_obj.export_as_plain_txt("-", "|")
	f_out = options[:path_to_output] + "-std.txt"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("txt_intron_phases") then
	output = gene_alignment_obj.export_as_plain_txt("-", nil)
	f_out = options[:path_to_output] + "-intron-phase.txt"
	write_verbosely_output_to_file(f_out, output)
end
		
if options[:output_format_list].include?("txt_only_introns") then 
	output = gene_alignment_obj.export_as_plain_txt(" ", "|")
	f_out = options[:path_to_output] + "-spaces.txt"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("txt_phylo") then 
	output = gene_alignment_obj.export_as_binary_alignment
	f_out = options[:path_to_output] + "-phylo.fas"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("svg") then 
	output = gene_alignment_obj.export_as_svg( options[:svg_options] )
	f_out = options[:path_to_output] + ".svg"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("pdb") then 
	# pdb option creates two scripts
	output1, output2 = gene_alignment_obj.export_as_pdb( options[:pdb_options] )
	f_out = options[:path_to_output] + "-color_exons.py"
	write_verbosely_output_to_file(f_out, output1)
	f_out = options[:path_to_output] + "-color_splicesites.py"
	write_verbosely_output_to_file(f_out, output2)
end

if options[:output_format_list].include?("stats") then 
	output = gene_alignment_obj.export_as_statistics(taxonomy_obj) # taxonomy_obj might be nil
	f_out = options[:path_to_output] + "-stats.txt"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("extensive_tax") then 
	output = gene_alignment_obj.export_as_taxonomy(taxonomy_obj)
	f_out = options[:path_to_output] + "-taxonomy.txt"
	write_verbosely_output_to_file(f_out, output)
end

puts " done."
