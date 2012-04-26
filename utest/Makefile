UTEST_CHECK = 1
TOPDIR	= ..
include $(TOPDIR)/Makefile.system

TARGET=openblas_utest

CUNIT_URL=http://downloads.sourceforge.net/project/cunit/CUnit/2.1-2/CUnit-2.1-2-src.tar.bz2
CUNIT_DIR=$(CURDIR)/CUnit-2.1-2

CUNIT_LIB=$(CUNIT_DIR)/lib/libcunit.a

CFLAGS+=-I$(CUNIT_DIR)/include

OBJS=main.o test_rot.o test_swap.o test_axpy.o test_dotu.o test_rotmg.o test_dsdot.o test_amax.o

all : run_test

CUnit-2.1-2-src.tar.bz2:
ifeq ($(OSNAME), Darwin)
	curl -O $(CUNIT_URL)
else
	wget $(CUNIT_URL)
endif

$(CUNIT_DIR): CUnit-2.1-2-src.tar.bz2
	@if test `$(MD5SUM) CUnit-2.1-2-src.tar.bz2 | $(AWK) '{print $$1}'` = 31c62bd7a65007737ba28b7aafc44d3a; then \
		echo $(TAR) xjf $< ;\
		$(TAR) xjf $< ; \
	else \
		rm -rf $(CUNIT_DIR) ;\
		echo "	Cannot download CUnit-2.1-2-src.tar.bz2 or the MD5 check sum is wrong (Please use orignal)."; \
		exit 1; \
	fi


$(CUNIT_LIB): $(CUNIT_DIR)
	(cd $(CUNIT_DIR); CC=$(CC) CFLAGS="$(CFLAGS)" ./configure --prefix=$(CUNIT_DIR))
	$(MAKE) -C $(CUNIT_DIR)
	$(MAKE) -C $(CUNIT_DIR) install

$(TARGET): $(CUNIT_LIB) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ ../$(LIBNAME) $(CUNIT_LIB) $(EXTRALIB)

run_test: $(TARGET)
	./$(TARGET)

clean:
	-rm -f *.o $(TARGET) 
	-rm -rf $(CUNIT_DIR)

libs:

