#!/usr/bin/perl

use strict;
use warnings;

use autodie;

use feature qw(say);
use File::Basename;

use Statistics::R;

use Data::Dumper;

#USAGE: perl create_heatmap_from_tsv.pl ~/git/PGMLab/data/cristina/

$|=1;

my $working_dir = $ARGV[0];
$working_dir =~ s/\/$//; # Remove a trailing slash
my $working_dirname = (split /\//, $working_dir)[-1];
my $output_tsv = "$working_dir/$working_dirname.pp.tsv";

open(my $tsv_fh, ">", $output_tsv);

opendir(my $dh, $working_dir);
my @sub_dirs = grep {-d "$working_dir/$_" && ! /^\.{1,2}$/} readdir($dh);

my $heatmap_filepath = "$working_dir/pp_heatmap.png";

my $header_flag = 0;
foreach my $sub_dir (@sub_dirs) {
    my $pp_filepath = "$working_dir/$sub_dir/$sub_dir.pp.tsv";
    next unless(-R $pp_filepath);
    my ($filename, $dir, $suffix) = fileparse($pp_filepath, qr"\..[^.]*$");
    my $pathway_name = (split /\./, $filename)[0];

    open(my $pp_fh, "<", $pp_filepath);
    
    my $header = <$pp_fh>;
    chomp $header;

    if ($header_flag == 0) {
        print $tsv_fh "$header\tpathway\n";
        $header_flag = 1;
    }
       
    my @column_names = split /\t/, $header;
    
    while(my $line = <$pp_fh>) {
        chomp $line;
        my @states = split /\t/, $line;
        my $gene = shift @states;
        
        #check to see if they are all equal
        unless (keys %{{ map {$_, 1} @states }} == 1) {
           print $tsv_fh "$line\t$pathway_name\n";
        }
    }
    close($pp_fh);
}  
   
close($tsv_fh);

my $R = Statistics::R->new();

$R->run(q'library(gplots)');
$R->run(qq'hmap<- read.table("$output_tsv", header = TRUE, sep = \'\t\', row.names = 1)');
$R->run(q'dim(hmap)');
$R->run(q'y<- as.matrix(hmap[1:8,1:344])');
$R->run(q'is.numeric(y)'); # it should be
$R->run(q'my_palette <- colorRampPalette(c("light green",  "purple"))(n = 199)');
$R->run(q'row_annotation <- (c (rep("yellow",2), rep("cyan", 1), rep("grey", 1), rep("black", 1), rep("blue", 2), rep("dark red", 1)))');
$R->run(q'heatmap.2(y, scale="none", key=TRUE,key.title=NA,keysize =1.00, key.xlab=NA, symkey=FALSE, density.info="none", trace="none", col=my_palette,cexRow=0.75, cexCol=0.35,margins =c(6,12), RowSideColors=c(row_annotation))');
$R->run(q'par(lend = 1)');
$R->run(q'legend("topright",      # location of the legend on the heatmap plot
          legend = c( "ERBB2", "Mitotic S","PIP3", "PTK6","RAF_MAP", "Rho_GTP"), # category labels
          col = c("yellow", "cyan", "grey", "black", "blue", "dark red"), border=FALSE, bty="n", lwd = 6, y.intersp = 0.7, cex=0.5)');

$R->run(q'dev.off()');
