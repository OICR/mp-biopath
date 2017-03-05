#!/usr/bin/perl

use strict;
use warnings;

use autodie;

use feature qw(say);
use Data::Dumper;

use LWP::UserAgent;  
use HTTP::Request;
use JSON;

## USAGE clear; perl get_expression_atlas_data.pl 2 gene-to-expression.tsv

# This will output the gene to expression values into the specified tsv file
# It will use the second "2" tissue type listed when running.

my $type_index = $ARGV[0] - 1;
my $out_file = $ARGV[1];

my $pathway_id_file = "../../PGM/list_of_pathways_to_be_exported.tsv";

open(my $pfh, "<", $pathway_id_file);

my @pathway_ids;
while (my $row = <$pfh>) {
   chomp $row;
   my @cells = split /\t/, $row;
   push @pathway_ids, $cells[0];
}

close($pfh);

my $URL_base = "https://www.ebi.ac.uk/gxa/widgets/heatmap/referenceExperiment?geneQuery=R-HSA-";

my @tissue_types;
my $count = 0;

my %gene_to_expression_value;

foreach my $pathway_id (@pathway_ids) {
	$count++;
	my $URL = "$URL_base$pathway_id";

	my $ua = LWP::UserAgent->new(ssl_opts => { verify_hostname => 1 });  
	my $header = HTTP::Request->new(GET => $URL);  
	my $request = HTTP::Request->new('GET', $URL, $header);  
	my $response_json = $ua->request($request);  

	if ($response_json->is_success){  
		my $data = decode_json($response_json->content);

		if ($count == 1) {
			my $headers = $data->{columnHeaders};
			for my $header (@$headers) {
				push @tissue_types, $header->{factorValue};
			}
		}

		my $rows = $data->{profiles}{rows};
		for my $row (@$rows) {
			my $name = $row->{name};
			my $expressions = $row->{expressions};
			my $node_expression_data = $expressions->[$type_index];
			if (defined $node_expression_data->{value}) {
				my $value = $node_expression_data->{value};
				$gene_to_expression_value{$name} = $value;
			}
			else {
				#say "no value";
			}
		}
	}
	elsif ($response_json->is_error){  
		print "Error:$URL\n";  
		print $response_json->error_as_HTML;  
	}
}

my $selected_tissue = $tissue_types[$type_index];
say "Selected Tissue Type:\t$selected_tissue";
print Dumper \@tissue_types;

open(my $fh, ">", $out_file); 

while(my($gene, $expression_value) = each %gene_to_expression_value) {
	say $fh "$gene\t$expression_value";
}
close($fh);
