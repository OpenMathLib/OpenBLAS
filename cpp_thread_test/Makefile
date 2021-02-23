include ../Makefile.rule

all :: dgemv_tester dgemm_tester omatcopy_rt_bench

dgemv_tester :
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemv_thread_safety.cpp ../libopenblas.a -lpthread -o dgemv_tester
	./dgemv_tester

dgemm_tester : dgemv_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemm_thread_safety.cpp ../libopenblas.a -lpthread -o dgemm_tester
	./dgemm_tester

omatcopy_rt_bench :
	$(CC) $(COMMON_OPT) -Wall -Wextra -Wshadow omatcopy_rt_bench.c -o omatcopy_rt_bench
	./omatcopy_rt_bench
clean ::
	rm -f dgemv_tester dgemm_tester omatcopy_rt_bench
