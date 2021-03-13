include ../Makefile.rule

all :: dgemv_tester dgemm_tester dpotrf_bug2562

dgemv_tester :
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemv_thread_safety.cpp ../libopenblas.a -lpthread -o dgemv_tester
	./dgemv_tester

dgemm_tester : dgemv_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow -fopenmp -std=c++11 dgemm_thread_safety.cpp ../libopenblas.a -lpthread -o dgemm_tester
	./dgemm_tester

dpotrf_bug2562 : dgemv_tester
	$(CXX) $(COMMON_OPT) -Wall -Wextra -Wshadow  -std=c++11 dpotrf_bug2562.cpp ../libopenblas.a -lpthread -o dpotrf_bug2562
	./dpotrf_bug2562

clean ::
	rm -f dgemv_tester dgemm_tester
