# !/usr/bin/env ruby

if ! RUBY_VERSION.include?("2.0") then
	Helper.abort( "You are using Ruby #{RUBY_VERSION}. For GenePainter, Ruby >= 2.0 is needed.")
end

# require gems and library
# TODO remove ruby-debug
require 'ruby-debug'
Dir[File.join(File.dirname(__FILE__), 'lib', '*.rb')].each {|file| require File.absolute_path(file) }

### parse command line arguments
args = ARGV.map {|ele| ele.dup}
OptParser.backwards_compability_check(args)
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

# make sure all _needed_ aligned sequences are of same length
common_aligned_seqs = Sequence.ensure_seqs_have_same_length(aligned_seqs, aligned_seqs_names, common_names)

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
		# TODO
		data_obj = GffToGene.new(data, gene_name)

	else
		Helper.abort("Cannot parse gene structure. Unknown file type #{file}.")
	end

	# create a gene object with this structure and sequence
	gene_obj = data_obj.to_gene # method to_gene returns a gene object containing the structure
	gene_obj.add_aligned_seq(common_aligned_seqs[ind]) # ... and aligned sequence
	gene_obj.add_alignment_range(options[:range] || []) # ... and range(s) of interesting parts [numbers match to the aligned sequence]
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
puts "Prepare output ... "
output_arr = []
f_out_extension = ""

### ... for every requested output format
options[:output_format].each do |format|
	case format
	when "alignment_with_intron_phases"
		output_arr = gene_alignment_obj.export_as_alignment_with_introns
		f_out_extension = ".fas"
		is_alignment = true
	when "txt_simple"
		GeneAlignment.exon_intron_placeholder=["-", "|"]
		output_arr = gene_alignment_obj.export_as_plain_txt
		f_out_extension = "-std.txt"
	when "txt_intron_phases"
		GeneAlignment.exon_intron_placeholder=["-", nil]
		output_arr = gene_alignment_obj.export_as_plain_txt
		f_out_extension = "-intron-phase.txt"
	when "txt_only_introns"
		GeneAlignment.exon_intron_placeholder=[" ", "|"]
		output_arr = gene_alignment_obj.export_as_plain_txt
		f_out_extension = "-spaces.txt"
	when "txt_phylo"
		GeneAlignment.exon_intron_placeholder=["0", "1"]
		output_arr = gene_alignment_obj.export_as_binary_alignment
		f_out_extension = "-phylo.fas"
		is_alignment = true
	when "svg"
		GeneAlignment.exon_intron_placeholder=["-", "|"]
		gene_alignment_obj.export_as_svg( options[:svg_options] )
		f_out_extension = ".svg"
	when "pdb"
		GeneAlignment.exon_intron_placeholder=["-", nil]
		output_arr = gene_alignment_obj.export_as_pdb( options[:pdb], options[:consensus] )
	else
		# this should never be executed, but it does not harm anyway
		puts "---"
		puts "Unknown output option. Provide plain text output instead."
		GeneAlignment.exon_intron_placeholder=["-", "|"]
		f_out_extension = "-std.txt"
		output_arr = gene_alignment_obj.export_as_plain_txt
	end
			
	print "\t"
	### add merged/consensus profile statistics 
	if format != "pdb" then
		if options[:merge] then
			print "calculating merged profile ... "
			stats_obj = GeneAlignment::Statistics.new(output_arr, is_alignment || false )
			merged_pattern = stats_obj.get_merged_exon_intron_pattern(is_alignment || false )
			output_arr << merged_pattern
		end
		if options[:consensus] then
			print "calculating consensus profile ... "
			stats_obj = GeneAlignment::Statistics.new(output_arr, is_alignment || false )
			consensus_pattern = stats_obj.get_consensus_exon_intron_pattern(options[:consensus], is_alignment || false )
			output_arr << consensus_pattern
		end
	end

	### output the output :-)
	f_out = options[:path_to_output] + f_out_extension
	print "writing output to #{f_out} ... "
	IO.write( f_out, output_arr.join("\n"), :mode => "w" )
	puts "done."

end



# TODO
# replace buggy Needleman-Wunsch algo by a working one
# inspect the code from marcel: /fab8/mahe/work/pyscipio/src/_fastnw (but there, one of the strings is DNA, the other protein and frameshifts are considered: in trace back, rucksprung-adresse is also stored)

