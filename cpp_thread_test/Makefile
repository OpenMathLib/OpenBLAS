include ../Makefile.rule

all :: dgemv_tester dgemm_tester test_gemm_omp

dgemv_tester :
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemv_thread_safety.cpp ../libopenblas.a -lpthread -o dgemv_tester
	./dgemv_tester

dgemm_tester : dgemv_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemm_thread_safety.cpp ../libopenblas.a -lpthread -o dgemm_tester
	./dgemm_tester

test_gemm_omp :
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp test_gemm_omp.cc ../libopenblas.a -lpthread -o test_gemm_omp
	./test_gemm_omp

clean ::
	rm -f dgemv_tester dgemm_tester test_gemm_omp
