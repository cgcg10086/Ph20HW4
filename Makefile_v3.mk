# comments 
# target: dependencies
#	(tab) action or no action 

# This is V3: it uses pattern rules and variables. 
# Need to imprve the python script to make this better. 

SRC_FILE = make_fig.py
REPORT_SRC=ph20_hw4_GeChen
OBJS=explicit implicit phasespace symplectic

.PHONY : figure1
	figure1: $(OBJS) $(SRC_FILE) 
		python $(SRC_FILE) $(OBJS) $@ 

.PHONY : figure2
	figure2: $(SRC_FILE) 
		python $(SRC_FILE) $(OBJS) $@ 

.PHONY : gitlog
	gitlog: 
		git log > git.log 

.PHONY : latex
	latex: 
		pdflatex $(REPORT_SRC) 

.PHONY : clean
	clean :
	        rm *.png *.aux *.fdb_latexmk *.fls *.log *.pdf 
