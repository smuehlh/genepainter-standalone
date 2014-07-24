# !/usr/bin/env ruby

require 'optparse'
require 'yaml'
require 'net/http'

require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'helper.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'sequence.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'taxonomyToGene.rb')

#default parameter  
Scipio_params = {:blattile => 7, :minid => 90, :maxmis => 7, :min_score => 0.3, :exhaust_align_size => 15000}
Webscipio_server = "www.webscipio.org"
Webscipio_port = "80"

# checks if species is searchable (= has genomes) via the webscipio api
# also returns the genomes found
def is_species_valid(species_query)
	url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches")
	post_parameters = {'search_genomes' => 'true', 'query' => species_query}
	response = Net::HTTP.post_form(url, post_parameters)
	id = response.body

	url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches/#{id}.yaml")
	response = Net::HTTP.get_response(url)
	yaml_string = response.body
	genomes = YAML::load(yaml_string)

	return genomes.any?, genomes
end

def select_best_genome(genomes)
	if genomes.size == 1 then
		# it's just one genome, simply use it! 
		return genomes.first
	else
		# select for: type, version, size
		# sort genomes by reverse order: worst genome is listed first, best last

		# in reverse order: from worst to best case
		genome_assembly_types = ["reads", "contigs", "supercontigs", "ultracontigs", "chromosome"]
		sorted_genomes = genomes.sort_by do |hash|
			ind_assembly_type =  genome_assembly_types.index(hash[:type]) || 0 # if its not found in type, assume worst type
			version = [ hash[:major_version] || 0,
				hash[:minor_version] || 0, 
				hash[:mini_version] || 0] # assume worst version if its not specified
			size = hash[:size] || 0 # assume worst size if its not speciified

			# sort by type, version, size
			[ind_assembly_type, version, size]
		end

		return sorted_genomes.last
	end
end

def run_scipio(genome, fasta)
	url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches")
	post_parameters = {'scipio_run' => 'true', 'target_file_path' => genome[:path], 'query' => fasta}
	post_parameters.merge!(Scipio_params)
	response = Net::HTTP.post_form(url, post_parameters)
	id = response.body

	url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches/#{id}.yaml")
	scipio_status_result = ["running", ""] # scipio returns "running" until its done!
	while( scipio_status_result[0] == "running" ) do 
		response = Net::HTTP.get_response(url)
		yaml_string = response.body
		scipio_status_result = YAML::load(yaml_string)
		sleep(5)
	end

	scipio_status = scipio_status_result[0]
	scipio_result = YAML::load(scipio_status_result[1])

	return scipio_status, scipio_result

end

def is_scipio_run_successful(scipio_status)
	return scipio_status == "finished"
end

def parse(args)

	# default options
	options = Hash.new

	opt_parser = OptionParser.new do |opts|

		opts.banner = "\nGenePainter v.2.0 maps gene structures onto multiple sequence alignments"
		opts.separator "Contact: Martin Kollmar (mako[at]nmr.mpibpc.mpg.de)"

		opts.separator ""
		opts.separator "Generate YAML-formatted gene structures for multiple sequence alignment (with WebScipio)"
		opts.separator "Usage: generate_yaml_for_multiple_species.rb -s <species_names_file> -i <fasta_file> [options]"
		opts.separator ""
		opts.separator "Gene structures are saved in files named by corresponding fasta headers"
		opts.separator ""

		opts.separator ""
		opts.on("-s", "--species-to-fasta <species_names_file>", String, 
			"Text-based file mapping fasta header to species names", 
			"Mandatory line-format:",
			"Fastaheader1[,Fastaheader2]:Species",
			"A list of all species available can be found at http://www.webscipio.org/webscipio/genome_list") do |file|
			Helper.file_exist_or_die(file)
			options[:path_to_species_mapping_file] = file
		end
		opts.on("-i", "--input <fasta_file>", String, 
			"Path to fasta-formatted multiple sequence alignment.") do |file|

			Helper.file_exist_or_die(file)
			headers, seqs = Sequence.read_in_alignment(file)
			if seqs.empty? then 
				Helper.abort "Invalid file #{file}. Must be fasta-formatted."
			end
			options[:fasta_sequences] = seqs
			options[:fasta_headers] = headers
		end
		opts.on("-p", "--path <path_to_output_files>", String, 
			"Path to folder where gene structures should be stored.") do |path|
			Helper.dir_exist_or_die(path)
			options[:output_path] = path
		end
		opts.on("-w", "--whitelist <path_to_whitelist>", String,
			"Text-based file containing a selection of fasta headers included in <species_names_file>", 
			"Only for fasta headers specified here gene structures are generated",
			"If not specified, gene structures are generated for all fasta headers listed in <species_names_file>",
			"Mandatory file-format:",
			"Each line contains exactly one fasta header exactly as it appers in <species_names_file>") do |file|
			Helper.file_exist_or_die(file)
			options[:fasta_headers_to_genenerate_strucutures_for] = IO.readlines(file, "\n").map{|line| line.chomp}
		end

		opts.separator ""
		opts.on_tail("-h", "--help", "Show this message") do 
			puts opts
			exit
		end

	end # optionparser

	if args.empty? then
		# display help and exit if program is called without any argument
		puts opt_parser.help
		exit
	end

	opt_parser.parse(args)

	# parse species to fasta file
	options[:species_with_corresponding_genes] = TaxonomyToGene.map_genenames_to_speciesnames(
		options[:path_to_species_mapping_file], options[:fasta_headers])
	options.delete(:path_to_species_mapping_file)


	# make sure mandatory arguments are present
	if options[:species_with_corresponding_genes].nil? || options[:fasta_sequences].empty?  then 
		Helper.abort "Missing mandatory argument. Please specify -s <species_names_file> and -i <fasta_file>." 
	end

	return options

	# use the own format of fatal error messages!				
	rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption => exc
		Helper.abort exc
end # parse()

args = ARGV
options = parse(args)
best_genome_by_species = {}

begin 

	## scipio-search in "best" genome for protein
	puts  "\nStart processing ... "
	generated_files = []

	options[:species_with_corresponding_genes].each do |species, genes|

		print species

		# select best genome
		species_api = species.gsub(" ", "_")
		valid_species, genomes = is_species_valid(species_api)
		if valid_species then
			## find "best" genome for species 
			best_genome = select_best_genome(genomes)
		end

		if best_genome.nil? || best_genome[:path].nil? || best_genome[:path].empty? then 
			Helper.warn "Cannot select a genome for species #{species}."
			print " ... no genome found!\n"
			next
		end

		# select genes
		genes_to_process = genes # default: generate gene structure for all genes
		if options[:fasta_headers_to_genenerate_strucutures_for] then 
			# a whitelist exists, only generate gene structures for those genes in whitelist
			genes_to_process = genes & options[:fasta_headers_to_genenerate_strucutures_for]
		end
		if genes_to_process.empty? then 
			# none of genes belonging to this species should be processed
			print " ... no gene structure to generate!\n"
			next
		else
			print "... generating gene structures for #{genes_to_process.join(", ")}"
		end

		# collect sequences mapped to this species
		fasta_this_species = []
		genes_to_process.each do |gene|
			if ind = options[:fasta_headers].index(gene) then 
				seq = options[:fasta_sequences][ind]
				fasta_this_species << Sequence.convert_strings_to_fasta(gene, seq)
			end
		end

		status, scipio_results = run_scipio(best_genome, fasta_this_species.join("\n"))

		if is_scipio_run_successful(status) then
			scipio_results.each do |fasta_header, result|
				f_name = "#{fasta_header}.yaml"
				if options[:output_path] then 
					f_name = File.join(options[:output_path], f_name)
				end
				fh = File.new(f_name, "w")
				fh.puts YAML::dump(result)
				fh.close

				generated_files.push f_name
			end
		else
			Helper.warn "#{fasta_header}: Nothing found with standard parameters. Visit www.webscipio.org to try different parameters."
		end 
	end

	puts " ... done."

	puts "Saved results to files #{generated_files.join(", ")}."

rescue => e
	Helper.abort "Something (#{e}) went wrong."
	exit 1
end
exit 0