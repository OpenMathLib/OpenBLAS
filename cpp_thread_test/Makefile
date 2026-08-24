TOPDIR = ..
include $(TOPDIR)/Makefile.system

all :: dgemv_tester dgemm_tester dgemm_mixed_tester

ifeq ($(OSNAME), WINNT)
all :: dgemm_shutdown_tester
endif

CALLBACK_THREAD_TEST =
ifeq ($(SMP),1)
ifeq ($(USE_OPENMP),1)
CALLBACK_THREAD_TEST = 1
else
ifneq ($(OSNAME),$(filter $(OSNAME),WINNT CYGWIN_NT))
CALLBACK_THREAD_TEST = 1
endif
endif
endif

ifeq ($(CALLBACK_THREAD_TEST),1)
.PHONY : dgemm_mixed_callback_tester
all :: dgemm_mixed_callback_tester
endif

dgemv_tester :
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -std=c++11 dgemv_thread_safety.cpp ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) -o dgemv_tester
	./dgemv_tester

dgemm_tester : dgemv_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -std=c++11 dgemm_thread_safety.cpp ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) -o dgemm_tester
	./dgemm_tester

dgemm_mixed_tester : dgemm_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -std=c++11 dgemm_thread_safety_mixed.cpp ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) -o dgemm_mixed_tester
	./dgemm_mixed_tester

dgemm_shutdown_tester : dgemm_mixed_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -std=c++11 dgemm_thread_safety_shutdown.cpp ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) -o dgemm_shutdown_tester
	./dgemm_shutdown_tester

ifeq ($(CALLBACK_THREAD_TEST),1)
dgemm_mixed_callback_tester : dgemm_mixed_tester
	./dgemm_mixed_tester --callback
endif

clean ::
	rm -f dgemv_tester dgemm_tester dgemm_mixed_tester dgemm_shutdown_tester
