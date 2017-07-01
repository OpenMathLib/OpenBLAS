TOPDIR  = ..
include $(TOPDIR)/Makefile.system



SRC = $(wildcard src/*.c)
OBJS = $(SRC:%.c=%.o)

TEST_SUITS = \
	slauum dlauum clauum zlauum \
	spotrf dpotrf cpotrf zpotrf \
	spbtrf dpbtrf cpbtrf zpbtrf \
	ssygst dsygst chegst zhegst \
	ssytrf dsytrf csytrf chetrf zsytrf zhetrf \
	sgetrf dgetrf cgetrf zgetrf \
	sgbtrf dgbtrf cgbtrf zgbtrf \
	strsyl dtrsyl ctrsyl ztrsyl \
	stgsyl dtgsyl ctgsyl ztgsyl \
	sgemmt dgemmt cgemmt zgemmt
TESTS = $(TEST_SUITS:%=test/%.pass)  # dummies
TEST_EXES = $(TEST_SUITS:%=test/%.x)

LINK_TEST = -L$(TOPDIR) -lopenblas -lgfortran -lm

.SECONDARY: $(TEST_EXES)
.PHONY: test

# ReLAPACK compilation

libs:	$(OBJS)
	@echo "Building ReLAPACK library $(LIBNAME)"
	$(AR) -rv  $(TOPDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(TOPDIR)/$(LIBNAME)

%.o: %.c config.h
	$(CC) $(CFLAGS) -c $< -o $@


# ReLAPACK testing

test: $(TEST_EXES) $(TESTS)
	@echo "passed all tests"

test/%.pass: test/%.x
	@echo -n $*:
	@./$< > /dev/null && echo " pass" || (echo " FAIL" && ./$<)

test/s%.x: test/x%.c test/util.o $(TOPDIR)/$(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=s $< test/util.o -o $@ $(LINK_TEST) $(TOPDIR)/$(LIBNAME) $(LINK_TEST)

test/d%.x: test/x%.c test/util.o $(TOPDIR)/$(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=d $< test/util.o -o $@ $(LINK_TEST) $(TOPDIR)/$(LIBNAME) $(LINK_TEST)

test/c%.x: test/x%.c test/util.o $(TOPDIR)/$(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=c $< test/util.o -o $@ $(LINK_TEST) $(TOPDIR)/$(LIBNAME) $(LINK_TEST)

test/z%.x: test/x%.c test/util.o $(TOPDIR)/$(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=z $< test/util.o -o $@ $(LINK_TEST) $(TOPDIR)/$(LIBNAME) $(LINK_TEST)


# cleaning up

clean:
	rm -f  $(OBJS) test/util.o test/*.x
