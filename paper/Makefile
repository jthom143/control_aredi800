# Generate pdf with pandoc


paper_outline.pdf : paper_outline.md
	pandoc paper_outline.md --latex-engine=xelatex --variable mainfont="Helvetica" --filter pandoc-fignos  -o $@