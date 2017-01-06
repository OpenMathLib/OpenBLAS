include ../make.inc

.SUFFIXES: .o .f
all: testlsame testslamch testdlamch testsecond testdsecnd testieee testversion

testlsame: lsame.o lsametst.o
	$(LOADER) $(LOADOPTS) -o $@ lsame.o lsametst.o

testslamch: slamch.o lsame.o slamchtst.o
	$(LOADER) $(LOADOPTS) -o $@ slamch.o lsame.o slamchtst.o

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(LOADER) $(LOADOPTS) -o $@ dlamch.o lsame.o dlamchtst.o

testsecond: second_$(TIMER).o secondtst.o
	@echo "[INFO] : TIMER value: $(TIMER) (given by make.inc)"
	$(LOADER) $(LOADOPTS) -o $@ second_$(TIMER).o secondtst.o

testdsecnd: dsecnd_$(TIMER).o dsecndtst.o
	@echo "[INFO] : TIMER value: $(TIMER) (given by make.inc)"
	$(LOADER) $(LOADOPTS) -o $@ dsecnd_$(TIMER).o dsecndtst.o

testieee: tstiee.o
	$(LOADER) $(LOADOPTS) -o $@ tstiee.o

testversion: ilaver.o LAPACK_version.o
	$(LOADER) $(LOADOPTS) -o $@ ilaver.o LAPACK_version.o

clean:
	rm -f *.o
.f.o:
	$(FORTRAN) $(OPTS) -c -o $@ $<

slamch.o: slamch.f ; $(FORTRAN) $(NOOPT) -c -o $@ $<
dlamch.o: dlamch.f ; $(FORTRAN) $(NOOPT) -c -o $@ $<
