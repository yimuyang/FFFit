
COP = gfortran 

#COP = f90 -C -g2 -r8 

PROG = fitff

#PROG = gibbs.dbg

DOTO = fitff.o        \
       numgrad.o       \
       analgrad.o         \
       deviation.o           \
       function.o          \
	   distance.o			\
	   energy.o

$(PROG):  $(DOTO)
	$(COP) $(DOTO) -o $(PROG) 


fitff.o: fitff.f90
	$(COP) -c fitff.f90	
	
distance.o: distance.f90
	$(COP) -c distance.f90

numgrad.o: numgrad.f90
	$(COP) -c numgrad.f90

analgrad.o: analgrad.f90
	$(COP) -c analgrad.f90

deviation.o: deviation.f90
	$(COP) -c deviation.f90

function.o: function.f90
	$(COP) -c function.f90

energy.o: energy.f90
	$(COP) -c energy.f90







