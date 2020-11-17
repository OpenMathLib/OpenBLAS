include ../Makefile.rule

all :: dgemv_tester dgemm_tester test_dgemm_omp

dgemv_tester :
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemv_thread_safety.cpp ../libopenblas.a -lpthread -o dgemv_tester
	./dgemv_tester

dgemm_tester : dgemv_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemm_thread_safety.cpp ../libopenblas.a -lpthread -o dgemm_tester
	./dgemm_tester

test_dgemm_omp:
	$(CC) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp test_dgemm_omp.c ../libopenblas.a -lpthread -o test_dgemm_omp
	./test_dgemm_omp
clean ::
	rm -f dgemv_tester dgemm_tester
