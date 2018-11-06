# comments 
# target: dependencies
#	(tab) action or no action 

# This is V1: the most basic version. There are lots of deplications. 


.PHONY : figure1
	figure1.png: 
		python make_fig.py argument1 argument2 

.PHONY : figure2
	figure2.png: 
		python make_fig.py argument3 argument4 

.PHONY : gitlog
	gitlog: 
		git log > git.log 

.PHONY : latex
	latex: 
		pdflatex hw4_GeChen.tex

.PHONY : clean
	clean :
	        rm *.png *.aux *.fdb_latexmk *.fls *.log *.pdf 
