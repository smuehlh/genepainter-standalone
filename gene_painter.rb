# !/usr/bin/env ruby

require 'yaml'

# require .rb files in library (including all subfolders)
Dir[File.join(File.dirname(__FILE__), 'lib', '**', '*.rb')].each do |file|
	require File.absolute_path(file)
end

# require ruby version 2.X
version_numbers = RUBY_VERSION.split(".")
major_version = version_numbers[0].to_i
if major_version < 2 then
        Helper.abort( "You are using Ruby #{RUBY_VERSION}. For GenePainter, Ruby >= 2.0 is needed.")
elsif major_version > 2
        Helper.warn("You are using Ruby #{RUBY_VERSION}. Please make sure that your Ruby is backwards compatible with Ruby 2.0, on which GenePainter was developed.")
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

print "Read in alignment ..."
## alignment
aligned_seqs_names, aligned_seqs = Sequence.read_in_alignment(options[:path_to_alignment])
puts  " done."

## genes
gene_names = [] # a list of all gene names
Dir.glob( File.join(options[:path_to_genestruct], "*.*") ) do |file|
	f_name = File.basename(file, ".*")
	gene_names << f_name
end

## initialize 
gene_objects = [] # a list of all gene objects
common_names = aligned_seqs_names & gene_names
genes_selected_for_output = [] # list of all gene names which are selected for output or data analysis
selection_criterium = nil

if common_names.empty? then 
	Helper.abort "No matches between alignment and gene structures."
end

## taxonomy
if options[:tax_options][:path_to_tax] && options[:tax_options][:path_to_tax_mapping] then
	## taxonomy, only if neccessary
	print "Read in taxonomy ..."
	taxonomy_obj = TaxonomyToGene.new( 
		options[:tax_options][:path_to_tax], 
		options[:tax_options][:path_to_tax_mapping], 
		options[:tax_options][:no_grep],
		options[:tax_options][:nice_grep],
		common_names )
	puts " done."

	if taxonomy_obj.species_not_found_in_NCBI_tax.any? then 
		puts "No taxonomy entry found for species #{taxonomy_obj.species_not_found_in_NCBI_tax.join(", ")}. Ignore them."
	end
	if taxonomy_obj.species_without_gene_structures.any? then 
		Helper.log "No gene structure for genes associated with species #{taxonomy_obj.species_without_gene_structures.join(", ")}. Ignore them."
	end
	if taxonomy_obj.get_all_species_linked_to_genes.empty? then 
		Helper.log "No taxonomy entry found for any species. Continue without taxonomy."
		taxonomy_obj = nil
	end
end

## create gene_objects, they will be in same order as they are in the alignment
print "Read in genes ..."

# handle selection criteria
if species = options[:select_by][:species] then 
	selection_criterium = "species #{species}"

	if taxonomy_obj && taxonomy_obj.species_with_corresponding_genes then
		# species to gene name mapping is already parsed, use it
		species2genes = taxonomy_obj.species_with_corresponding_genes
	else 
		# parse species-to-genename-mapping
		species2genes = TaxonomyToGene.map_genenames_to_speciesnames(options[:tax_options][:path_to_tax_mapping])
	end

	if species2genes[species] then 
		genes_selected_for_output = species2genes[species]
	else
		Helper.abort "Selection criterium species #{species} not found in file #{options[:tax_options][:path_to_tax_mapping]}"
	end
end
if list = options[:select_by][:list] then
	selection_criterium = "list #{list}" 
	genes_selected_for_output = list & common_names
	if (list - common_names).any? then 
		missing = list - common_names
		selection_criterium = "#{selection_criterium}.\nCould not find #{missing.join(", ")} amoung gene names."
	end
end
if regex = options[:select_by][:regex] then 
	selection_criterium = "regular expression #{regex.inspect}"
end
if selection_criterium.nil? then 
	# no selection criterum was specified, this means to use all data and output all data
	genes_selected_for_output = common_names
end

# make sure all _needed_ aligned sequences are of same length
common_aligned_seqs = Sequence.ensure_seqs_have_same_length(aligned_seqs, aligned_seqs_names, common_names)

# remove common gaps if neccessary
if options[:ignore_common_gaps] then 
	common_aligned_seqs = Sequence.remove_common_gaps(common_aligned_seqs, 
		{delete_all_common_gaps: true} # remove all common gaps
	) 
end


# special case: gap before/after intron
# correct intron position (default: before gap) if there is any other sequence with and intron at same pos as gap-end
introns_before_gap_pos_gene = {}

# merge gene structure and sequence (and taxonomy, if present) into a gene object
common_names.each_with_index do |gene_name, ind|

	# if selection is based on regexp, apply it to gene name
	if options[:select_by][:regex] then 
		genes_selected_for_output.push gene_name if options[:select_by][:regex].match(gene_name)
	end

	# analyse selection of all data only? then read in only selected genes!
	if options[:selection][:analyse_sel_output_sel] then 
		# skip all genes which do not belong to selection
		next if ! genes_selected_for_output.include?(gene_name)
	end

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
	gene_obj.add_taxonomy(taxonomy_obj.to_gene(gene_name)) if taxonomy_obj

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
	print "Correcting intron positions flanking alignment gaps ..."

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
# inform which data are used for output and analysis
if selection_criterium then 
	selection_type = "analyse selection and output selection" if options[:selection][:analyse_sel_output_sel]
	selection_type = "analyse all data and output selection" if options[:selection][:analyse_all_output_sel]
	selection_type = "analyse selection on basis of all data and output selection" if options[:selection][:analyse_sel_on_all_output_sel] 
	Helper.print_selection(genes_selected_for_output, selection_type, selection_criterium)
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

# removing all genes which should not be part of output
# all other genes are completely deleted from gene_alignment_obj
if options[:selection][:analyse_all_output_sel] || options[:selection][:analyse_sel_on_all_output_sel] then 
	if options[:selection][:analyse_all_output_sel] then 
		is_delete_introns_not_occuring_in_sel = false
	else 
		# options[:selection][:analyse_sel_on_all_output_sel] is true
		is_delete_introns_not_occuring_in_sel = true
	end
	
	print "Removing genes from output ..."
	if genes_selected_for_output.empty? then 
		Helper.abort "Selection criterium too restrictive. No gene selected for output."
	end
	gene_alignment_obj.reduce_gene_set_for_output(genes_selected_for_output,
		is_delete_introns_not_occuring_in_sel,
		options[:consensus], 
		options[:merge], 
		{ genes_belonging_to_selected_taxa: genes_belonging_to_selected_taxa || [], 
			is_intron_exclusive_for_selected_taxa: options[:tax_options][:is_exclusive] || false 
		}
	)
	puts "done."
end

### prepare output for every requested format
# checking for each possible format, if it was specified, is faster than looping through all specified ones
# NEVER change names of output files! (for webserver)
puts "Prepare output ... "
if options[:output_format_list].include?("alignment_with_intron_phases") then 
	# this is in most cases the master format
	output = gene_alignment_obj.export_as_alignment_with_introns
	f_out = options[:path_to_output] + "-alignment.fas"
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

if options[:output_format_list].include?("txt_fuzzy") then 
	output = gene_alignment_obj.export_as_plain_text_with_fuzzy_intron_pos("-", "|", options[:fuzzy_window])
	f_out = options[:path_to_output] + "-fuzzy.txt"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("svg") then 
	# svg creates normal, reduced or both output files
	output1, output2 = gene_alignment_obj.export_as_svg( options[:svg_options] )
	if output1 then 
		f_out = options[:path_to_output] + "-normal.svg"
		write_verbosely_output_to_file(f_out, output1)
	end
	if output2 then 
		f_out = options[:path_to_output] + "-reduced.svg"
		write_verbosely_output_to_file(f_out, output2)	
	end
	if options[:svg_options][:generate_merged_pattern] then 
		# only needed for webserver, that uses a extra file containing only merged pattern
		output1, output2 = gene_alignment_obj.export_as_svg_only_merged_pattern( options[:svg_options] )
		if output1 then 
			f_out = options[:path_to_output] + "-normal-merged.svg"
			write_verbosely_output_to_file(f_out, output1)
		end
		if output2 then 
			f_out = options[:path_to_output] + "-reduced-merged.svg"
			write_verbosely_output_to_file(f_out, output2)
		end
	end
end

if options[:output_format_list].include?("pdb") then 
	# pdb option creates two scripts
	# no need to feed in gene names for which output should be generated
	output1, output2 = gene_alignment_obj.export_as_pdb( options[:pdb_options] )
	f_out = options[:path_to_output] + "-color_exons.py"
	write_verbosely_output_to_file(f_out, output1)
	f_out = options[:path_to_output] + "-color_splicesites.py"
	write_verbosely_output_to_file(f_out, output2)
end

if options[:output_format_list].include?("stats") then 
	output = gene_alignment_obj.export_as_statistics("-", "|")
	f_out = options[:path_to_output] + "-stats.txt"
	write_verbosely_output_to_file(f_out, output)
end

if options[:output_format_list].include?("extensive_tax") then 
	f_out = options[:path_to_output] + "-taxonomy.txt"
	is_success = catch(:no_taxonomy) do
		output = gene_alignment_obj.export_as_taxonomy
		write_verbosely_output_to_file(f_out, output)
		true
	end
	if ! is_success then 
		puts "\t taxonomy is missing: cannot write to #{f_out}"
	end
end

if options[:tax_options][:generate_list_intron_positios_per_taxon] then 
	# only needed for webserver
	f_out = options[:path_to_output] + "-taxonomy-intron-numbers.txt"
	is_success = catch(:no_taxonomy) do 
		output = gene_alignment_obj.export_as_taxonomy_list_of_intron_positions_per_taxon_only
		write_verbosely_output_to_file(f_out, output)
		true
	end
	if ! is_success then 
		puts "\t taxonomy is missing: cannot write to #{f_out}"
	end
end

if options[:output_format_list].include?("tree") then 
	f_out_phb = options[:path_to_output] + "-tree.phb"
	is_success = catch(:no_taxonomy) do 
		output = gene_alignment_obj.export_as_tree
		write_verbosely_output_to_file(f_out_phb, output)
		f_out_svg = options[:path_to_output] + "-tree.svg"
		path_to_python_script = File.join(File.dirname(__FILE__), 'tools', 'phb2svg.py')
		is_success_python = system 'python', *[path_to_python_script, '--spacing=5', '--minstep=80', '--textwidth=150', f_out_phb, f_out_svg]
		if ! is_success_python then 
			Helper.log "Could not convert #{f_out_phb} to SVG."
		else
			puts "\t writing output to #{f_out_svg} ... "
		end
		true
	end
	if ! is_success then 
		puts "\t taxonomy is missing: cannot write to #{f_out_phb}"
	end
end

puts " done."
