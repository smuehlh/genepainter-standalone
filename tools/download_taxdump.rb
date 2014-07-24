# !/usr/bin/env ruby

require 'optparse'
require 'timeout'
require File.join(File.expand_path("..",File.dirname(__FILE__)), 'lib', 'helper.rb')

def parse(args)

	options = Hash.new
	options[:name] = "taxdump"

	opt_parser = OptionParser.new do |opts|

		opts.banner = "\nGenePainter v.2.0 maps gene structures onto multiple sequence alignments"
		opts.separator "Contact: Martin Kollmar (mako[at]nmr.mpibpc.mpg.de)"

		opts.separator ""
		opts.separator "Download NCBI taxonomy dump to gene_painter directory"
		opts.separator "Usage: download_taxdump.rb [-o <file_name>]"
		opts.separator ""

		opts.separator ""
		opts.separator "Options:"
		opts.on("-o", "--outfile <file_name>", String, 
			"Name of the output file") do |file|
			options[:name] = File.basename(file, ".*")
		end

		opts.separator ""
		opts.on_tail("-h", "--help", "Show this message") do 
			puts opts
			exit
		end

	end # optionparser

	opt_parser.parse(args)

	return options

	# use the own format of fatal error messages!				
	rescue OptionParser::MissingArgument, OptionParser::InvalidArgument, OptionParser::InvalidOption => exc
		Helper.abort exc
end # parse()

def call_wget(f_path)
	wget_path = `which wget`.chomp
	system(wget_path, "-O", f_path, "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
end

args = ARGV
options = parse(args)

max_secs = 60*10 # wait max. 10 minutes
output_path = File.join(File.expand_path('../..',__FILE__), options[:name] + ".tar.gz") # download to parent dir of this file 
begin
        status = Timeout::timeout(max_secs) { call_wget(output_path) }
rescue Timeout::Error => exc
        Helper.abort(exc)
end
if ! status then 
	$stderr.puts "\nFatal error: wget failed"
	exit 1
end
puts "Downloaded NCBI taxonomy to file #{output_path}"
exit 0