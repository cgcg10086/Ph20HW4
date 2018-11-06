# comments 
# target: dependencies
#	(tab) action or no action 

# This is V2: it uses automatic variables $@, $^, $< to aovid repeatings. 


.PHONY : figure1
	figure1.png: 
		python make_fig.py argument1 $@ 

.PHONY : figure2
	figure2.png: 
		python make_fig.py argument3 $@ 

.PHONY : gitlog
	gitlog: 
		git log > git.log 

.PHONY : latex
	latex: 
		pdflatex hw4_GeChen.tex

.PHONY : clean
	clean :
	        rm *.png *.aux *.fdb_latexmk *.fls *.log *.pdf 
