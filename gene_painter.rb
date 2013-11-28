# !/usr/bin/env ruby

if ! RUBY_VERSION.include?("2.0") then
	Helper.abort( "You are using Ruby #{RUBY_VERSION}. For GenePainter, Ruby >= 2.0 is needed.")
end

# require gems and library
# TODO remove ruby-debug
require 'ruby-debug'
Dir[File.join(File.dirname(__FILE__), 'lib', '*.rb')].each {|file| require File.absolute_path(file) }

args = ARGV.map {|ele| ele.dup}
options = OptParser.parse(args)

### open log file and specify its automatic closure at_exit
fh_log = File.open(options[:path_to_log], "w")
at_exit { fh_log.close }
fh_log.puts "GenePainter logfile"
fh_log.puts "Called with arguments: #{ARGV.join(" ")}"
fh_log.puts "Date: #{Time.now}"

### read in data, use only intersect of alingment and gene structures
print "Read in alignment and genes ..."
## alignment
aligned_seqs_names, aligned_seqs = [], []
IO.foreach(options[:path_to_alignment]) do |line|
	line.chomp!
	if line[0] == ">" then
		# fasta header
		aligned_seqs_names << line[1..-1] 
	else
		# fasta sequence

		n_seqs = aligned_seqs_names.size
		# a new sequence or another line for the last sequence?
		if aligned_seqs.size < n_seqs then
			# new sequence
			aligned_seqs << line
		else
			# add to last sequence
			aligned_seqs[n_seqs - 1] += line # -1: convert number of elements n_seqs to an index
		end
	end
end
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
common_names.each do |gene_name|

	# find gene structure file with belongs to this gene
	matching_files = Dir.glob( File.join(options[:path_to_genestruct], "#{gene_name}.*") ) 

	# double-check that the gene was matched un-ambiguously to gene structure files
	if matching_files.size > 1 then
		Helper.abort("Cannot parse gene structure. #{gene_name} is ambiguous file name.")
	else
		file = matching_files[0]
		f_extension = File.extname(file)
		aligned_seq_ind = aligned_seqs_names.index(gene_name)
	end

	# read in file and parse gene structure
	data_obj = nil
	data = IO.read(file)
	if f_extension.casecmp(".yaml") == 0 then
		# yaml format
		data_obj = YamlToGene.new(data, gene_name)

	elsif f_extension.casecmp(".gff") == 0 then
		# gff format
		# TODO
		data_obj = GffToGene.new(data, gene_name)

	else
		Helper.abort("Cannot parse gene structure. Unknown file type #{file}.")
	end

	# create a gene object with this structure
	gene_obj = data_obj.to_gene # method to_gene returns a gene object
	gene_obj.aligned_seq = aligned_seqs[aligned_seq_ind] # ... and aligned sequence

	gene_objects << gene_obj
end

puts " done."

# inform which data are not used
Helper.print_intersect_and_diff_between_alignment_and_gene(aligned_seqs_names, gene_names, fh_log)

### align genes
puts ""
print "Aligning genes ..."
gene_alignment_obj = GeneAlignment.new(gene_objects)
puts " done."

### prepare output 
print "Prepare output ... "
output_str = ""
f_out_extension = ""

case options[:output_format]
when "alignment_with_intron_phases"
	output_str = gene_alignment_obj.export_as_alignment_with_introns
	f_out_extension = ".fas"

when "txt_simple", "txt_intron_phases", "txt_only_introns"

	# how to display exons and introns?
	# possible combinations:
	# txt_simple: exon="-",intron="|"
	# txt_intron_phases: exon="-",intron=["0"|"1"|"2"]
	# txt_only_introns: exon=" ",intron="|"
	exon_representation_in_output = "-" # gap symbol
	intron_representation_in_output = nil # intron phases 
	case options[:output_format]
	when "txt_simple"
		intron_representation_in_output = "|"
	when "txt_only_introns"
		exon_representation_in_output = " "
		intron_representation_in_output = "|"
	end
	f_out_extension = ".txt"
	output_str = gene_alignment_obj.export_as_plain_txt(exon_representation_in_output, intron_representation_in_output)

when "txt_phylo"
	output_str = gene_alignment_obj.export_as_binary_alignment
	f_out_extension = ".fas"
else
	puts ""
	puts "*** No!"
end
		
### output the output :-)
print "writing output ... "
IO.write( options[:path_to_output] + f_out_extension, output_str, :mode => "w" )
puts "done."


