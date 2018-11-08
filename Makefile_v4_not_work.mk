# comments 
# target: dependencies
#	(tab) action or no action 

# This is V3: it uses pattern rules and variables. 
# Need to imprve the python script to make this better. 

SRC_FILE = make_fig.py
REPORT_SRC=ph20_hw4_GeChen
OBJS=explicit implicit phasespace symplectic

all: report # defalut target 

.PHONY : report
	report: $(OBJS) # report depands on the figures (for which numerical integration methods in hw3)
		git log > git.log 
		pdflatex $(REPORT_SRC) 

.PHONY : $(OBJS) 
	$(OBJS): $(SRC_FILE) # figures depands on the makefig.py script 
		python $< $@ 

.PHONY : clean # delete all generated files in the directory 
	clean :
	        rm *.png *.aux *.fdb_latexmk *.fls *.log *.pdf 
