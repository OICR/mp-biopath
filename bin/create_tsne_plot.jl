using TSne  #to install run command Pkg.clone("git://github.com/lejon/TSne.jl.git")
using DataFrames

# This is from https://github.com/lejon/TSne.jl/blob/35bad1c5cfd78d44e5d89bbca89a85555943b9db/examples/demo-csv.jl

# USAGE: julia bin/create_tsne_plot.jl haveheader --labelcolname=pathway_label /home/awright/git/PGM/gecco/mpbiopath_results.summary.for_tsne.tsv

doc = """Use t-SNE to generate a PDF called myplot.pdf from an input CSV file. Default assumption is to have no header and no labels. If these are available in the CSV these must be given as arguments.
Usage:
  demo-csv.jl <filename>
  demo-csv.jl [--labelcol=<col>] <filename>
  demo-csv.jl haveheader [--labelcol=<col>] <filename>
  demo-csv.jl [--labelcolname=<colname>] <filename>
  demo-csv.jl haveheader [--labelcolname=<colname>] <filename>
Options:
  -h --help     Show this screen.
  --version     Show version.
  --filename=   Path to CSV file
  --noheader    The CSV file does not have a header row
  --nolabel    The CSV file does not have a label column
"""

function normalize(A)
	for col in 1:size(A)[2]
        	std(A[:,col]) == 0 && continue
        	A[:,col] = (A[:,col]-mean(A[:,col])) / std(A[:,col])
	end
	A
end

using DocOpt

arguments = docopt(doc, version=v"0.0.1")
dump(arguments)

if nothing==arguments["--labelcol"]
	lblcol = -1
else
	lblcol = parse(Int64,arguments["--labelcol"])
end

df = readtable(arguments["<filename>"], header = nothing!=arguments["haveheader"], separator = '\t')

if nothing!=arguments["--labelcolname"]
	lblcol = find(x -> x==symbol(arguments["--labelcolname"]),names(df))[1]
end

println("Data is $df")
if lblcol>0
	labels = df[:,lblcol]
end

dataset = df[filter(x -> x!=lblcol,1:ncol(df)),]
data = float(convert(Array,dataset))
# Normalize the data, this should be done if there are large scale differences in the dataset
X = normalize(data)

# Run t-SNE
Y = tsne(X, 2, 50, 1000, 20.0)

using Gadfly
if lblcol>0
	theplot = plot(x=Y[:,1], y=Y[:,2], color=labels)
else
	theplot = plot(x=Y[:,1], y=Y[:,2])
end
draw(PDF("tsne-plot.pdf", 8inch, 6inch), theplot)
