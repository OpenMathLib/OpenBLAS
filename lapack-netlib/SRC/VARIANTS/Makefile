include ../../make.inc

#######################################################################
#  This is the makefile to create a the variants libraries for LAPACK.
#  The files are organized as follows:
#       CHOLRL -- Right looking block version of the algorithm, calling Level 3 BLAS
#       CHOLTOP -- Top looking block version of the algorithm, calling Level 3 BLAS
#       LUCR -- Crout Level 3 BLAS version of LU factorization
#       LULL -- left-looking Level 3 BLAS version of LU factorization
#       QRLL -- left-looking Level 3 BLAS version of QR factorization
#       LUREC -- an iterative version of Sivan Toledo's recursive LU algorithm[1].  
#       For square matrices, this iterative versions should
#       be within a factor of two of the optimum number of memory transfers.
#
# [1] Toledo, S. 1997. Locality of Reference in LU Decomposition with
#  Partial Pivoting. SIAM J. Matrix Anal. Appl. 18, 4 (Oct. 1997),
#  1065-1081. http://dx.doi.org/10.1137/S0895479896297744
#######################################################################

VARIANTSDIR=LIB

CHOLRL = cholesky/RL/cpotrf.o cholesky/RL/dpotrf.o cholesky/RL/spotrf.o cholesky/RL/zpotrf.o

CHOLTOP = cholesky/TOP/cpotrf.o cholesky/TOP/dpotrf.o cholesky/TOP/spotrf.o cholesky/TOP/zpotrf.o

LUCR = lu/CR/cgetrf.o lu/CR/dgetrf.o lu/CR/sgetrf.o lu/CR/zgetrf.o

LULL = lu/LL/cgetrf.o lu/LL/dgetrf.o lu/LL/sgetrf.o lu/LL/zgetrf.o

LUREC = lu/REC/cgetrf.o lu/REC/dgetrf.o lu/REC/sgetrf.o lu/REC/zgetrf.o

QRLL = qr/LL/cgeqrf.o qr/LL/dgeqrf.o qr/LL/sgeqrf.o qr/LL/zgeqrf.o  qr/LL/sceil.o


all: cholrl choltop lucr lull lurec qrll

cholrl: $(CHOLRL)
	$(ARCH) $(ARCHFLAGS) $(VARIANTSDIR)/cholrl.a $(CHOLRL)
	$(RANLIB) $(VARIANTSDIR)/cholrl.a

choltop: $(CHOLTOP)
	$(ARCH) $(ARCHFLAGS) $(VARIANTSDIR)/choltop.a $(CHOLTOP)
	$(RANLIB) $(VARIANTSDIR)/choltop.a

lucr: $(LUCR)
	$(ARCH) $(ARCHFLAGS) $(VARIANTSDIR)/lucr.a $(LUCR)
	$(RANLIB) $(VARIANTSDIR)/lucr.a

lull: $(LULL)
	$(ARCH) $(ARCHFLAGS) $(VARIANTSDIR)/lull.a $(LULL)
	$(RANLIB) $(VARIANTSDIR)/lull.a

lurec: $(LUREC)
	$(ARCH) $(ARCHFLAGS) $(VARIANTSDIR)/lurec.a $(LUREC)
	$(RANLIB) $(VARIANTSDIR)/lurec.a
	
qrll: $(QRLL)
	$(ARCH) $(ARCHFLAGS) $(VARIANTSDIR)/qrll.a  $(QRLL)
	$(RANLIB) $(VARIANTSDIR)/qrll.a


.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
	
clean:
	rm -f $(CHOLRL) $(CHOLTOP) $(LUCR) $(LULL) $(LUREC) $(QRLL) \
	      $(VARIANTSDIR)/*.a