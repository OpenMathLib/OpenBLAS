#
# tests of installed headers and libs
#

INSTALLDIR = ../install

.PHONY: all
all: test_sched_include

test_sched_include: test_sched_include.c
	$(CC) -c -I$(INSTALLDIR)/include $<

.PHONY: clean
clean:
	rm -f *.o
