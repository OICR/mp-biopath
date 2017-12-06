#!/usr/bin/env perl

use strict;
use warnings;

use autodie;
use feature qw(say);

open(my $in_fh, "<", "./data/essential_9606_all.txt");

<$in_fh>;  #header

my %gene_status;
while(my $line = <$in_fh>) {
   chomp($line);

   my @fields = split /\t/, $line;

   if ($fields[5] eq "TextMining") {
       next;
   }

   $gene_status{$fields[2]}{$fields[3]} = 1;
}

close($in_fh);

open(my $out_fh, ">", "./data/essential_9606_all_gene_status.txt");

foreach my $gene (keys %gene_status) {
    next unless (exists $gene_status{$gene}{E});

    my $essentiality = (exists $gene_status{$gene}{"NE"})? "Conditional": "Essential";
    print $out_fh "$gene\t$essentiality\n";
}
