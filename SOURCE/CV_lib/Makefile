
# include directories list
# IDIR = /home/amit/Softwares/GSL1.15/gsl-1.15/include
IDIR = /usr/gapps/qmc/libs/INTEL/gsl-1.9/include 

CC = icc 

AR = ar

ARFLAGS = rvu 

RANLIB = ranlib

CFLAGS = -g -O2 -I$(IDIR)
#CFLAGS = -g -O0 -C -I$(IDIR)

LIBS = -lgsl  -lgslcblas -lm 

# LFLAGS = -L/home/amit/Softwares/GSL1.15/gsl-1.15/lib/
LFLAGS = #-L/usr/gapps/qmc/libs/INTEL/gsl-1.9/lib

OBJ = main.o 

LIBOBJ = CV.o 

MYEXEC = main

.PHONY: clean c

lib : libcv.a

$(MYEXEC) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(MYEXEC) $(LFLAGS) $(LIBS)

$(OBJ) : main.c
	$(CC) $(CFLAGS) -c main.c

libcv.a   : $(LIBOBJ)
	rm -f libcv.a
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

c:clean

clean:
	rm -f *.o *.a *~ \#* main 
# DO NOT DELETE

CV.o: CV.c Common.c Nlist.c Ql.c
