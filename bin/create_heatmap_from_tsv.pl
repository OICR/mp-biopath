use strict;
use warnings;

use autodie;
use feature qw(say);

use Statistics::R;
use Data::Dumper;

#USAGE: perl bin/create_heatmap_from_tsv.pl ~/git/PGM/gecco/mpbiopath-gecco.results.tsv ~/git/PGM/list_of_pathways_to_be_exported.tsv

my $data_file = $ARGV[0];
my $pathway_file = $ARGV[1];

open(my $pl_fh, "<", $pathway_file);

my %pathway_color_map;
while (my $line = <$pl_fh>) {
    chomp $line;
    my @parts = split /\t/, $line;
    if (defined $parts[3]) {
        $pathway_color_map{$parts[2]} = $parts[3];
    }
}

my $output_file = "heatmap";

open(my $fh, "<", $ARGV[0]);

my $header = <$fh>;
chomp $header;
my @column_names = split /\t/, $header;

my %labels;
my $column_index = 0;
foreach my $column_name (@column_names) {
    $labels{$1}{index} = $column_index if ($column_name =~ /(.*)_label/);
    $column_index++;
}

my $number_of_obs = 0;
my %label_values;
my @label_names = keys %labels;
my %pathways;
my @color_column;
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split /\t/, $line;
    foreach my $label (@label_names) {
        my $label_value = $fields[$labels{$label}{index}];
        if (exists $labels{$label}{values}) {
            push @{$labels{$label}{values}}, $label_value;
        }
        else {
            $labels{$label}{values} = [$label_value];
        }
        my $color = $pathway_color_map{$label_value};
        push @color_column, "\"$color\"";
        $pathways{"\"$label_value\""} = "\"$color\"";
    }
    $number_of_obs++;
}

close($fh);

my $cases = scalar(@column_names) - scalar(@label_names) - 1;

my @pathway_list = keys %pathways;
my $pathway_string = join ", ", @pathway_list;

my @color_list = values %pathways; ##need to get vlaues for the pwathway
my $color_string = join ", ", @color_list;

my $color_column_string = join ", ", @color_column;

my $R = Statistics::R->new();
$R->run(q'library(gplots)');
$R->run(qq'hmap<- read.table("$data_file", header = TRUE, sep = "\t", row.names = 1)');
$R->run(qq'y <- as.matrix(hmap[1:$number_of_obs,1:$cases])');
$R->run(q'y <- y[] + 0.01');
$R->run(q'my_palette <- colorRampPalette(c("blue", "white",  "red"))(n = 199)');
$R->run(qq'row_annotation <- (c($color_column_string))');
$R->run(q'par(lend = 1)');
$R->run(qq'svg("./$output_file.svg")');
$R->run(q'heatmap.2(log10(y),
                    scale="none",
                    key=TRUE,
                    key.title=NA,
                    keysize =1.00,
                    key.xlab=NA,
                    symkey=FALSE,
                    density.info="none",
                    col=my_palette,
                    trace="none",
                    cexRow=0.3,
                    cexCol=0.05,
                    margins=c(12,15),
                    RowSideColors=c(row_annotation))');
$R->run(qq'legend("topright",      # location of the legend on the heatmap plot
           legend = c($pathway_string), # category labels
           col = c($color_string),
           border=FALSE, bty="n", lwd = 6, y.intersp = 0.7, cex=0.5)');
$R->run(qq'svg("$output_file.svg")');
$R->run(q'dev.off()');
