UTEST_CHECK = 1
TOPDIR	= ..
include $(TOPDIR)/Makefile.system

TARGET=openblas_utest
CUNIT_LIB=/usr/local/lib/libcunit.a

OBJS=main.o test_rot.o

all : run_test

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ ../$(LIBNAME) $(CUNIT_LIB)

run_test: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o $(TARGET)

libs:

