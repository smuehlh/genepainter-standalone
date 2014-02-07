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
	if output.end_with?("\n") then 
		output << "\n"
	end
	IO.write( f_name, output, :mode => "w" )
end

### open log file, and specify its automatic closure at_exit
Helper.open_log( options[:path_to_log] )
at_exit { Helper.close_log }
Helper.log "Program call: #{ARGV.join(" ")}"

### read in data, use only intersect of alignment and gene structures
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
common_aligned_seqs = Sequence.remove_common_gaps(common_aligned_seqs, {is_alignment: true}) if options[:ignore_common_gaps]

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
#	gene_obj.add_alignment_range(options[:range] || []) # ... and range(s) of interesting parts [numbers match to the aligned sequence]
	gene_obj.reduce_gene_to_range(options[:range]) if options[:range].any?
	gene_objects << gene_obj
end

puts " done."

# inform which data are not used
Helper.print_intersect_and_diff_between_alignment_and_gene(aligned_seqs_names, gene_names)

# select subset of genes belonging to certain taxa
if options[:tax_options].any? then 
	puts ""
	print "Preparing taxonomy ... "
	taxonomy_obj = Taxonomy.new(options[:tax_options][:path_to_tax], options[:tax_options][:path_to_tax_mapping])
	genes_within_taxa = taxonomy_obj.get_genenames_belonging_to_selected_taxa(common_names, options[:tax_options][:selected_taxa])
	puts " done."
	Helper.print_genes_within_taxa(genes_within_taxa, options[:tax_options][:selected_taxa])
end

### align genes
puts ""
print "Aligning genes ..."

# initiate an gene alignment object 
# this calculates kind of an "master format", the exon-intron-patterns for each gene
gene_alignment_obj = GeneAlignment.new(gene_objects, 
	options[:consensus], 
	options[:merge], 
	{ genes_within_taxa: genes_within_taxa, is_exclusive: options[:tax_options][:is_exclusive] }
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


puts " done."
