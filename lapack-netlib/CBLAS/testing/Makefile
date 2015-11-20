#
# The Makefile compiles c wrappers and testers for CBLAS.  
#

include ../../make.inc

# Archive files necessary to compile
LIB = $(CBLASLIB) $(BLASLIB) 

# Object files for single real precision
stestl1o = c_sblas1.o

stestl2o = c_sblas2.o c_s2chke.o auxiliary.o c_xerbla.o

stestl3o = c_sblas3.o c_s3chke.o auxiliary.o c_xerbla.o

# Object files for double real precision
dtestl1o = c_dblas1.o

dtestl2o = c_dblas2.o c_d2chke.o auxiliary.o c_xerbla.o

dtestl3o = c_dblas3.o c_d3chke.o auxiliary.o c_xerbla.o

# Object files for single complex precision
ctestl1o = c_cblas1.o

ctestl2o = c_cblas2.o c_c2chke.o auxiliary.o c_xerbla.o

ctestl3o = c_cblas3.o c_c3chke.o auxiliary.o c_xerbla.o

# Object files for double complex precision
ztestl1o = c_zblas1.o 

ztestl2o = c_zblas2.o c_z2chke.o auxiliary.o c_xerbla.o

ztestl3o = c_zblas3.o c_z3chke.o auxiliary.o c_xerbla.o

all: all1 all2 all3 
all1: stest1 dtest1 ctest1 ztest1 
all2: stest2 dtest2 ctest2 ztest2
all3: stest3 dtest3 ctest3 ztest3

clean:
	rm -f core *.o a.out x*
cleanobj:
	rm -f core *.o a.out 
cleanexe:
	rm -f x* 

stest1: xscblat1
dtest1: xdcblat1
ctest1: xccblat1
ztest1: xzcblat1

stest2: xscblat2
dtest2: xdcblat2
ctest2: xccblat2
ztest2: xzcblat2

stest3: xscblat3
dtest3: xdcblat3
ctest3: xccblat3
ztest3: xzcblat3

#
# Compile each precision
#

# Single real
xscblat1: $(stestl1o) c_sblat1.o
	$(LOADER) $(LOADOPTS) -o xscblat1 c_sblat1.o $(stestl1o) $(LIB)
xscblat2: $(stestl2o) c_sblat2.o
	$(LOADER) $(LOADOPTS) -o xscblat2 c_sblat2.o $(stestl2o) $(LIB)
xscblat3: $(stestl3o) c_sblat3.o
	$(LOADER) $(LOADOPTS) -o xscblat3 c_sblat3.o $(stestl3o) $(LIB)
# Double real
xdcblat1: $(dtestl1o) c_dblat1.o
	$(LOADER) $(LOADOPTS) -o xdcblat1 c_dblat1.o $(dtestl1o) $(LIB)
xdcblat2: $(dtestl2o) c_dblat2.o
	$(LOADER) $(LOADOPTS) -o xdcblat2 c_dblat2.o $(dtestl2o) $(LIB)
xdcblat3: $(dtestl3o) c_dblat3.o
	$(LOADER) $(LOADOPTS) -o xdcblat3 c_dblat3.o $(dtestl3o) $(LIB)
 
# Single complex
xccblat1: $(ctestl1o) c_cblat1.o
	$(LOADER) $(LOADOPTS) -o xccblat1 c_cblat1.o $(ctestl1o) $(LIB)
xccblat2: $(ctestl2o) c_cblat2.o
	$(LOADER) $(LOADOPTS) -o xccblat2 c_cblat2.o $(ctestl2o) $(LIB)
xccblat3: $(ctestl3o) c_cblat3.o
	$(LOADER) $(LOADOPTS) -o xccblat3 c_cblat3.o $(ctestl3o) $(LIB)

# Double complex 
xzcblat1: $(ztestl1o) c_zblat1.o
	$(LOADER) $(LOADOPTS) -o xzcblat1 c_zblat1.o $(ztestl1o) $(LIB)
xzcblat2: $(ztestl2o) c_zblat2.o
	$(LOADER) $(LOADOPTS) -o xzcblat2 c_zblat2.o $(ztestl2o) $(LIB)
xzcblat3: $(ztestl3o) c_zblat3.o
	$(LOADER) $(LOADOPTS) -o xzcblat3 c_zblat3.o $(ztestl3o) $(LIB)
   

# RUN TESTS
run:
	@echo "--> TESTING CBLAS 1 - SINGLE PRECISION <--"
	@./xscblat1 > stest1.out
	@echo "--> TESTING CBLAS 1 - DOUBLE PRECISION <--"
	@./xdcblat1  > dtest1.out
	@echo "--> TESTING CBLAS 1 - COMPLEX PRECISION <--"
	@./xccblat1  > ctest1.out
	@echo "--> TESTING CBLAS 1 - DOUBLE COMPLEX PRECISION <--"
	@./xzcblat1	 > ztest1.out
	@echo "--> TESTING CBLAS 2 - SINGLE PRECISION <--"
	@./xscblat2 < sin2 > stest2.out
	@echo "--> TESTING CBLAS 2 - DOUBLE PRECISION <--"
	@./xdcblat2 < din2 > dtest2.out
	@echo "--> TESTING CBLAS 2 - COMPLEX PRECISION <--"
	@./xccblat2 < cin2 > ctest2.out
	@echo "--> TESTING CBLAS 2 - DOUBLE COMPLEX PRECISION <--"
	@./xzcblat2	< zin2	> ztest2.out
	@echo "--> TESTING CBLAS 3 - SINGLE PRECISION <--"
	@./xscblat3 < sin3 > stest3.out
	@echo "--> TESTING CBLAS 3 - DOUBLE PRECISION <--"
	@./xdcblat3 < din3 > dtest3.out
	@echo "--> TESTING CBLAS 3 - COMPLEX PRECISION <--"
	@./xccblat3 < cin3 > ctest3.out
	@echo "--> TESTING CBLAS 3 - DOUBLE COMPLEX PRECISION <--"
	@./xzcblat3	< zin3 > ztest3.out	
	
.SUFFIXES: .o .f .c

.c.o:
	$(CC) -c $(CFLAGS) -I ../include -o $@ $<
	
.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
