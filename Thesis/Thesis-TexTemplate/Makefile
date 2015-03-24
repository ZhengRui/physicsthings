dvi: thesis.tex
	latex thesis
pdf: thesis.dvi
	gbk2uni thesis
	latex thesis
	dvipdfmx thesis
clean:
	find . -name '*.log' -exec rm -f {} \;  -o -name '*.aux' -exec rm -f {} \;
cleanall:
	rm -f thesis.aux thesis.dvi thesis.idx thesis.lof thesis.log thesis.lot thesis.out  thesis.thm thesis.toc
