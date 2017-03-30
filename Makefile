KEYOUTPUT_FILE=./data/keyoutputs.tsv

.PHONY: install
install: install-julia-packages

.PHONY: install-julia-packages
install-julia-packages:
	julia bin/INSTALL.jl
