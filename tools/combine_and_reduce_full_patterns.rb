# !/usr/bin/env ruby

require 'optparse'
require 'ruby-debug' # FIXME
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'sequence.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'helper.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'geneAlignment.rb')

# read in several, full exon-intron-patterns and combines them into single, reduced exon-intron pattern

def parse(args)

	options = Hash.new
	options[:input_files] = []
	options[:outfile] = nil
	options[:is_special_patterns_only] = false
	options[:special_pattern_names] = []

	opt_parser = OptionParser.new do |opts|

		opts.banner = "\nGenePainter v.2.0 maps gene structures onto multiple sequence alignments"
		opts.separator "Contact: Martin Kollmar (mako[at]nmr.mpibpc.mpg.de)"

		opts.separator ""
		opts.separator "Combine multiple, full exon-intron-patterns into single, reduced one."
		opts.separator "Usage: combine_and_reduce_full_patterns.rb -i <list_of_input_files> -o <file_name>"
		opts.separator ""

		opts.separator "Input: List or regex (including the path)"
		opts.on("-i", "--input x,y,z", Array, 
			"Combine patterns in files x,y,z", 
			"List must consist of at least one file") do |list|
			list = list.each {|file| Helper.file_exist_or_die(file) }
			options[:input_files] = list
		end
		opts.on("--input-regex <regex>", String,
			"OR specify regex for matching files, including the path to the files") do |opt|
			path = File.dirname(opt)
			file_regex = File.basename(opt)
			puts "Searching in dir #{path} for files matching #{file_regex} ..."
			list = []
			Dir.glob( File.join(path, "*") ).each do |file|
				if file.match(file_regex) then 
					list << file
				end
			end
			options[:input_files] = list.sort
		end

		opts.separator ""
		opts.separator "Output"
		opts.on("-o", "--outfile <file_name>", String, 
			"Name of the output file") do |file|
			options[:outfile] = file
		end

		opts.separator ""
		opts.separator "Options:"
		opts.on("--only-special-patterns x,y,z", Array, 
			"Combine special patterns only.", 
			"x,y,z may be 'Consensus' and 'Merged'") do |opt|
			options[:is_special_patterns_only] = true 

			if opt.include?("Consensus") then 
				options[:special_pattern_names].push "Consensus"
			elsif opt.include?("Merged") then 
				options[:special_pattern_names].push "Merged"
			end
		end

		opts.separator ""
		opts.on_tail("-h", "--help", "Show this message") do 
			puts opts
			exit
		end

	end # optionparser

	opt_parser.parse(args)

	if args.empty? then
		# display help and exit if program is called without any argument
		puts opt_parser.help
		exit
	end

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
@long_name_lenght = GeneAlignment.max_length_gene_name * 2
@is_use_long_names = false
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
	max_length_fname = (GeneAlignment.max_length_gene_name / 1.5).floor
	fname = File.basename(file)[0...max_length_fname]
	name = name.sub(">", "") if name.start_with?(">") # delete starting ">" in name
	chop_pattern_name("#{fname}_#{name}")
end

def generate_new_long_name(name, file)
	fname = File.basename(file)
	name = name.sub(">", "") if name.start_with?(">") # delete starting ">" in name
	@is_use_long_names = true
	"#{fname}_#{name}".ljust(@long_name_lenght)
end
def collect_exon_intron_patterns(file, is_special_patterns_only, special_pattern_names)
	arr = []

	IO.foreach(file) do |line|
		line = line.chomp
		name, pattern = split_in_name_and_pattern(line)
		# uniq_name = generate_new_name(name, file)
		uniq_name = generate_new_long_name(name, file)
		new_name_with_pattern = [uniq_name, pattern].join()

# only to trim phylo-input
# if ! pattern.include?("|") then 
# 	next
# end
		is_name_matches_special_pattern_names = false 
		special_pattern_names.each do |sp_name|
			if name.match(sp_name) then 
				is_name_matches_special_pattern_names = true 
			end
		end
		if !is_special_patterns_only || ( is_special_patterns_only && is_name_matches_special_pattern_names ) then 
			# use all patterns or only the special patterns (and in that case, this pattern is one of the special ones)
			arr.push new_name_with_pattern
		end
	end
	return arr
end

options = parse(ARGV)
path_to_log = File.join( File.dirname(options[:outfile]), "#{File.basename(options[:outfile], ".*")}.log" )
Helper.open_log( path_to_log )
at_exit { Helper.close_log }
Helper.log "Program call: #{ARGV.join(" ")}"

all_exon_intron_patterns = []
options[:input_files].each do |file|
 	all_exon_intron_patterns = all_exon_intron_patterns.concat collect_exon_intron_patterns(
 		file, options[:is_special_patterns_only], options[:special_pattern_names]
 		)
end
if @is_use_long_names then 
	start_col = @long_name_lenght
else
	start_col = GeneAlignment.max_length_gene_name
end

# options to trim phylo-input: delete_all_common_gaps: true, intron_pattern: false
reduced_patterns = Sequence.remove_common_gaps( all_exon_intron_patterns, {start_col: start_col, intron_pattern: true} )

IO.write( options[:outfile], reduced_patterns.join("\n"), :mode => "w" )
