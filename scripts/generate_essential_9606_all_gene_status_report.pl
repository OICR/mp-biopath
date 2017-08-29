#!/usr/bin/env perl

use strict;
use warnings;

use autodie;
use feature qw(say);

open(my $fh, "<", "./data/essential_9606_all_gene_status.txt");


my $number_essential = 0;
my $number_conditional = 0;
while(my $line = <$fh>) {
    chomp $line;

    my ($gene, $status) = split /\t/, $line;

    if ($status eq "Conditional") {
        $number_conditional++;
    } 
    else {
        $number_essential++;
    }
}
my $total = $number_conditional + $number_essential;

say "Conditional\t$number_conditional";
say "Essential\t$number_essential";
say "Total\t\t$total";

my $percent_essential = $number_essential/$total * 100;

printf("Percent Essential\t%2d%%\n", $percent_essential);

close($fh);
