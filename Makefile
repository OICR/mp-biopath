KEYOUTPUT_FILE=./data/keyoutputs.tsv

.PHONY: install
install: install-julia-packages

.PHONY: install-julia-packages
install-julia-packages:
	julia bin/INSTALL.jl

.PHONY: update-keyoutput-table
update-keyoutput-table:
	curl "https://docs.google.com/spreadsheets/d/1HgLlibdVW_qR2xvwcMT1v2QwipBf7BAoKw17cU0Wxx4/export?gid=0&format=tsv" > ${KEYOUTPUT_FILE}


