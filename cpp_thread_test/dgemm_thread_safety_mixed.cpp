#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>
#ifdef OPENBLAS_USE_GENERATED_CBLAS_H
#include "generated/cblas.h"
#else
#include "../cblas.h"
#endif
#include "cpp_thread_safety_common.h"

std::atomic<uint32_t> callbackInvocations(0);

void thread_callback(int sync, openblas_dojob_callback doJob, int numJobs,
                     size_t jobDataElementSize, void* jobData, int doJobData){
	(void)sync;
	callbackInvocations.fetch_add(1, std::memory_order_relaxed);
	std::vector<std::thread> workers;
	workers.reserve(numJobs);
	char* jobs = static_cast<char*>(jobData);
	for(int i=0; i<numJobs; i++)
		workers.emplace_back(doJob, i, jobs + i * jobDataElementSize, doJobData);
	for(auto& worker : workers)
		worker.join();
}

void compute_dgemm_pair(std::vector<double>& transA, std::vector<double>& noTransA, std::vector<double>& B, double* firstOutput, double* secondOutput, const blasint randomMatSize, const bool sameVariant){
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, randomMatSize, 2, 2, 1.0, &transA[0], randomMatSize, &B[0], 2, 0.0, firstOutput, 2);
	if (sameVariant)
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, randomMatSize, 2, 4, 1.0, &transA[0], randomMatSize, &B[0], 2, 0.0, secondOutput, 2);
	else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, randomMatSize, 2, 4, 1.0, &noTransA[0], 4, &B[0], 2, 0.0, secondOutput, 2);
}

void run_worker(std::vector<double>& transA, std::vector<double>& noTransA, std::vector<double>& B, const std::vector<double>& referenceFirst, const std::vector<double>& referenceSecond, const blasint randomMatSize, const uint32_t numTestRounds, const bool sameVariant, std::atomic<uint32_t>& readyThreads, std::atomic<bool>& startThreads, uint32_t& mismatches){
	std::vector<double> firstOutput(static_cast<size_t>(randomMatSize) * 2);
	std::vector<double> secondOutput(static_cast<size_t>(randomMatSize) * 2);
	const size_t outputBytes = static_cast<size_t>(randomMatSize) * 2 * sizeof(double);
	uint32_t localMismatches = 0;

	readyThreads.fetch_add(1);
	while (!startThreads.load())
		std::this_thread::yield();

	for(uint32_t i=0; i<numTestRounds; i++){
		compute_dgemm_pair(transA, noTransA, B, &firstOutput[0], &secondOutput[0], randomMatSize, sameVariant);
		if (std::memcmp(&firstOutput[0], &referenceFirst[0], outputBytes) != 0 ||
		    std::memcmp(&secondOutput[0], &referenceSecond[0], outputBytes) != 0)
			localMismatches++;
	}

	mismatches = localMismatches;
}

int main(int argc, char* argv[]){
	blasint randomMatSize = 262144;
	uint32_t numConcurrentThreads = 8;
	uint32_t numTestRounds = 200;
	uint32_t maxHwThreads = GetMaxHwThreads();
	bool sameVariant = false;
	bool useCallback = false;

	if (maxHwThreads < numConcurrentThreads)
		numConcurrentThreads = maxHwThreads;

	std::vector<std::string> positionalArgs;
	for (int i = 1; i < argc; i++){
		std::cout<<argv[i]<<std::endl;
		if (std::string(argv[i]) == "--callback")
			useCallback = true;
		else
			positionalArgs.push_back(argv[i]);
	}

	if (!positionalArgs.empty() && positionalArgs.size() != 3 && positionalArgs.size() != 4){
		std::cout<<"ERROR: expected: [<M> <threads> <rounds> [sameVariant]] [--callback]"<<std::endl;
		return 1;
	}

	if(!positionalArgs.empty()){
		randomMatSize = std::stoul(positionalArgs[0]);
		numConcurrentThreads = std::stoul(positionalArgs[1]);
		numTestRounds = std::stoul(positionalArgs[2]);
		if (positionalArgs.size() == 4)
			sameVariant = std::stoul(positionalArgs[3]) != 0;
	}

	FailIfThreadsAreZero(numConcurrentThreads);

	const size_t matrixElements = static_cast<size_t>(randomMatSize) * 4;
	const size_t outputElements = static_cast<size_t>(randomMatSize) * 2;
	std::vector<double> transA(matrixElements);
	std::vector<double> noTransA(matrixElements);
	std::vector<double> B(8);
	std::vector<double> referenceFirst(outputElements);
	std::vector<double> referenceSecond(outputElements);
	std::vector<std::thread> threads(numConcurrentThreads);
	std::vector<uint32_t> mismatchBlock(numConcurrentThreads);
	std::atomic<uint32_t> readyThreads(0);
	std::atomic<bool> startThreads(false);

	std::cout<<"*----------------------------------*\n";
	std::cout<<"| Mixed DGEMM thread safety tester |\n";
	std::cout<<"*----------------------------------*\n";
	std::cout<<"Tall-skinny DGEMM M dimension: "<<randomMatSize<<'\n';
	std::cout<<"Number of concurrent calls into OpenBLAS : "<<numConcurrentThreads<<'\n';
	std::cout<<"Number of testing rounds : "<<numTestRounds<<'\n';
	std::cout<<"Second DGEMM uses "<<(sameVariant ? "the same transpose variant" : "a different transpose variant")<<'\n';
	std::cout<<"OpenBLAS internal threads : "<<openblas_get_num_threads()<<'\n';
	if (useCallback)
		std::cout<<"Thread execution backend : caller callback\n";
	std::cout<<"This test will need "<<(static_cast<uint64_t>(matrixElements) * 2 * 8 + static_cast<uint64_t>(outputElements) * (2 + 2 * numConcurrentThreads) * 8)/static_cast<double>(1024*1024)<<" MiB of RAM\n"<<std::endl;

	std::cout<<"Filling matrices with deterministic values..."<<std::flush;
	for(size_t i=0; i<matrixElements; i++){
		transA[i] = static_cast<double>(i % 512);
		noTransA[i] = static_cast<double>((i * 7) % 512);
	}
	std::cout<<"done\n";
	std::cout<<"Filling RHS matrix..."<<std::flush;
	for(uint32_t i=0; i<8; i++){
		B[i] = 0.1 * static_cast<double>(i + 1);
	}
	std::cout<<"done\n";

	std::cout<<"Computing reference results..."<<std::flush;
	compute_dgemm_pair(transA, noTransA, B, &referenceFirst[0], &referenceSecond[0], randomMatSize, sameVariant);
	std::cout<<"done\n";

	if (useCallback)
		openblas_set_threads_callback_function(thread_callback);

	std::cout<<"Testing mixed CBLAS DGEMM thread safety\n";
	std::cout<<"Launching "<<numConcurrentThreads<<" worker threads..."<<std::flush;
	for(uint32_t i=0; i<numConcurrentThreads; i++){
		threads[i] = std::thread(run_worker, std::ref(transA), std::ref(noTransA), std::ref(B), std::cref(referenceFirst), std::cref(referenceSecond), randomMatSize, numTestRounds, sameVariant, std::ref(readyThreads), std::ref(startThreads), std::ref(mismatchBlock[i]));
	}
	while (readyThreads.load() != numConcurrentThreads)
		std::this_thread::yield();
	startThreads.store(true);
	std::cout<<"done\n";

	std::cout<<"Waiting for worker threads to finish..."<<std::flush;
	uint32_t mismatches = 0;
	for(uint32_t i=0; i<numConcurrentThreads; i++){
		threads[i].join();
		mismatches += mismatchBlock[i];
	}
	std::cout<<"done\n";

	if (useCallback) {
		const uint32_t invocations = callbackInvocations.load();
		std::cout<<"Thread callback invocations: "<<invocations<<std::endl;
		if (invocations == 0) {
			std::cout<<"Thread callback was not invoked!"<<std::endl;
			return 1;
		}
	}

	std::cout<<"Mixed DGEMM mismatches: "<<mismatches<<std::endl;
	if (mismatches != 0) {
		std::cout<<"Mixed CBLAS DGEMM thread safety test FAILED!"<<std::endl;
		return 1;
	}

	std::cout<<"Mixed CBLAS DGEMM thread safety test PASSED!\n"<<std::endl;
	return 0;
}
