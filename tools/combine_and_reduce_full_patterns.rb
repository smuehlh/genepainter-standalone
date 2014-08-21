# !/usr/bin/env ruby

require 'optparse'
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'sequence.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'helper.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'geneAlignment.rb')

# read in several, full exon-intron-patterns and combines them into single, reduced exon-intron pattern

def parse(args)

	options = Hash.new
	options[:input_files] = []
	options[:outfile] = nil
	options[:is_special_patterns_only] = false

	opt_parser = OptionParser.new do |opts|

		opts.banner = "\nGenePainter v.2.0 maps gene structures onto multiple sequence alignments"
		opts.separator "Contact: Martin Kollmar (mako[at]nmr.mpibpc.mpg.de)"

		opts.separator ""
		opts.separator "Combine multiple, full exon-intron-patterns into single, reduced one."
		opts.separator "Usage: combine_and_reduce_full_patterns.rb -i <list_of_input_files> -o <file_name>"
		opts.separator ""

		opts.separator ""
		opts.on("-i", "--input")

		opts.on("-i", "--input x,y,z", Array, 
			"Combine patterns in files x,y,z", 
			"List must consist of at least one file") do |list|
			list = list.each {|file| Helper.file_exist_or_die(file) }
			options[:input_files] = list
		end

		opts.on("-o", "--outfile <file_name>", String, 
			"Name of the output file") do |file|
			options[:outfile] = file
		end

		opts.separator "Options:"
		opts.on("--only-special-patterns", 
			"Combine merged and consensus pattern only") do |opt|
			options[:is_special_patterns_only] = opt
		end

		opts.separator ""
		opts.on_tail("-h", "--help", "Show this message") do 
			puts opts
			exit
		end

	end # optionparser

	opt_parser.parse(args)

	if options[:input_files].empty? then 
		Helper.abort "Invalid usage of '--input': Must contain at least one file"
	end
	if ! options[:outfile] then 
		Helper.abort "Invalid usage of '--outfile': Specify file name for output file"
	end

	return options

	# use the own format of fatal error messages!				
	rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption => exc
		Helper.abort exc
end # parse()

def get_prefix_of_additional_patterns
	prefixes = []
	prefixes.push GeneAlignment.merged_structure_name
	prefixes.push GeneAlignment.consensus_structure_name.sub(/_.*$/, "")
	prefixes.push GeneAlignment.taxonomy_structure_name
	return prefixes
end
def split_in_name_and_pattern(line)
	name, pattern = "", ""
	name = line[0...GeneAlignment.max_length_gene_name]
	pattern = line[GeneAlignment.max_length_gene_name..-1]
	return name, pattern
end
def chop_pattern_name(str)
	# add ">" at beginning of string
	# left justify string to certain length
	(">" << str)[0...GeneAlignment.max_length_gene_name].ljust(GeneAlignment.max_length_gene_name)
end
def generate_new_name(name, file)
	max_length_fname = (GeneAlignment.max_length_gene_name / 2.0).floor
	fname = File.basename(file)[0...max_length_fname]
	name = name.sub(">", "") if name.start_with?(">") # delete starting ">" in name
	chop_pattern_name("#{fname}_#{name}")
end
def collect_exon_intron_patterns(file, is_special_patterns_only)
	arr = []
	special_pattern_names = get_prefix_of_additional_patterns

	IO.foreach(file) do |line|
		line = line.chomp
		name, pattern = split_in_name_and_pattern(line)
		uniq_name = generate_new_name(name, file)
		new_name_with_pattern = [uniq_name, pattern].join()
		if !is_special_patterns_only || ( is_special_patterns_only && special_pattern_names.include?(name) ) then 
			# use all patterns or this pattern is one of the special patterns
			arr.push new_name_with_pattern
		end
	end
	return arr
end

options = parse(ARGV)

all_exon_intron_patterns = []
options[:input_files].each do |file|
 	all_exon_intron_patterns = all_exon_intron_patterns.concat collect_exon_intron_patterns(file, options[:is_special_patterns_only])
end

reduced_patterns = Sequence.remove_common_gaps( all_exon_intron_patterns, {start_col: GeneAlignment.max_length_gene_name} )

IO.write( options[:outfile], reduced_patterns.join("\n"), :mode => "w" )
