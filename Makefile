CPP = g++
FLAGS = -o
LFLAGS = -larmadillo

INDIR = ./input

vibrational_analysis: vibrational_analysis.cpp
	$(CPP) $(FLAGS) vibrational_analysis vibrational_analysis.cpp $(LFLAGS)

run_all: vibrational_analysis
	./vibrational_analysis $(INDIR)/CH4.txt
	./vibrational_analysis $(INDIR)/CO2.txt
	./vibrational_analysis $(INDIR)/H2.txt
	./vibrational_analysis $(INDIR)/H2O.txt
	./vibrational_analysis $(INDIR)/HF.txt
	./vibrational_analysis $(INDIR)/HO.txt
	./vibrational_analysis $(INDIR)/NH3.txt
	./vibrational_analysis $(INDIR)/O2.txt

clean:
	-rm ./vibrational_analysis