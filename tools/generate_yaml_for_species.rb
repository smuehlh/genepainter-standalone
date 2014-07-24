# !/usr/bin/env ruby

require 'optparse'
require 'yaml'
require 'net/http'
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'helper.rb')
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'sequence.rb')

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
	specified_yaml_name = nil # a temporary variable

	opt_parser = OptionParser.new do |opts|

		opts.banner = "\nGenePainter v.2.0 maps gene structures onto multiple sequence alignments"
		opts.separator "Contact: Martin Kollmar (mako[at]nmr.mpibpc.mpg.de)"

		opts.separator ""
		opts.separator "Generate YAML-formatted gene structure for one or more protein sequences of a given species (with WebScipio)"
		opts.separator "Usage: generate_yaml_for_species.rb -s <species_name> -i <fasta_file> [-p <path_to_output_files>]"
		opts.separator ""
		opts.separator "Gene structures are saved in files named by corresponding fasta headers"
		opts.separator ""

		opts.separator ""
		opts.on("-s", "--species \"<species_name>\"", String, 
			"Species encoding the specified protein", 
			"Species should be wrapped with \"\" to preserve spaces.",
			"A list of all species available can be found at http://www.webscipio.org/webscipio/genome_list") do |str|
			str_api = str.gsub(" ", "_")
			valid_species, genomes = is_species_valid(str_api)
			if ! valid_species then 
				Helper.abort( "No genome file found for species #{str}") 
			end
			options[:species] = str_api
			options[:available_genomes] = genomes
		end
		opts.on("-i", "--input <fasta_file>", String, 
			"Path to fasta-formatted protein sequence(s). Might be a multiple sequence alignment of sequences encoded by same species.") do |file|

			Helper.file_exist_or_die(file)
			headers, seqs = Sequence.read_in_alignment(file)
			if seqs.empty? then 
				Helper.abort "Invalid file #{file}. Must be fasta-formatted."
			end
			options[:fasta] = IO.read(file)
		end
		opts.on("-p", "--path <path_to_output_files>", String, 
			"Path to folder where gene structures should be stored.") do |path|
			Helper.dir_exist_or_die(path)
			options[:output_path] = path
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

	# make sure mandatory arguments are present
	if options[:species].nil? || options[:fasta].empty?  then 
		Helper.abort "Missing mandatory argument. Please specify -s <species_name> and -i <fasta_file>." 
	end

	# options.delete(:fasta_header) # not needed any more

	return options

	# use the own format of fatal error messages!				
	rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption => exc
		Helper.abort exc
end # parse()

args = ARGV
options = parse(args)

begin 

	## find "best" genome for species 
	best_genome = select_best_genome(options[:available_genomes])
	if best_genome.nil? || best_genome[:path].nil? || best_genome[:path].empty? then 
		Helper.abort "Cannot select a genome for species #{options[species]}."
	end
	puts "Selected genome #{File.basename(best_genome[:path])}."

	## scipio-search in "best" genome for protein

	print "Starting Scipio run ... "
	status, scipio_results = run_scipio(best_genome, options[:fasta])
	puts " done."

	if is_scipio_run_successful(status) then
		print "Saving results to file(s) "

		scipio_results.each do |fasta_header, result|
			f_name = "#{fasta_header}.yaml"
			if options[:output_path] then 
				f_name = File.join(options[:output_path], f_name)
			end
			fh = File.new(f_name, "w")
			fh.puts YAML::dump(result)
			fh.close

			print "#{f_name} "
		end

		puts "."
	else
		Helper.warn "#{fasta_header}: Nothing found with standard parameters. Visit www.webscipio.org to try different parameters."
	end 


rescue => e
	Helper.abort "Something (#{e}) went wrong."
	exit 1
end
exit 0