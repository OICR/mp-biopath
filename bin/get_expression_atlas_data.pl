#!/usr/bin/perl

use strict;
use warnings;

use autodie;

use feature qw(say);
use Data::Dumper;

use LWP::UserAgent;  
use HTTP::Request;
use JSON;

## USAGE clear; perl get_expression_atlas_data.pl gene-to-expression.tsv

# This will output the gene to expression values into the specified tsv file
# It will use the second "2" tissue type listed when running.

my $out_file = $ARGV[0];

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

my %tissue_to_expression_value;
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
	    my $i = 0;
	    foreach my $tissue_type (@tissue_types) {
                my $node_expression_data = $expressions->[$i];
                if (defined $node_expression_data->{value}) {
                     my $value = $node_expression_data->{value};
                     $tissue_to_expression_value{$i}{$name} = $value;
                }
                else {
                    #say "no value";
                }
		$i++;
            }
        }
    }
    elsif ($response_json->is_error){  
        print "Error:$URL\n";  
        print $response_json->error_as_HTML;  
    }
}

open(my $fh, ">", $out_file); 

while (my($tissue_type, $tissue_to_expression) = each %tissue_to_expression_value) {
    while (my($gene, $expression_value) = each %$tissue_to_expression) {
        say $fh "$tissue_type\t$gene\t$expression_value";
    }
}
close($fh);

open($fh, ">", "$out_file.index");

for my $i (0 .. $#tissue_types) {
    say $fh "$i\t$tissue_types[$i]";
}

close($fh);
