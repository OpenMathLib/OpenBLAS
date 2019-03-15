#include <iostream>
#include <vector>
#include <random>
#include <future>
#include <omp.h>
#include "/opt/OpenBLAS_zen_serial/include/cblas.h"

const blasint randomMatSize = 1024; //dimension of the random square matrices used
const uint32_t numConcurrentThreads = 52; //number of concurrent calls of the functions being tested
const uint32_t numTestRounds = 8; //number of testing rounds before success exit

inline void pauser(){
    /// a portable way to pause a program
    std::string dummy;
    std::cout << "Press enter to continue...";
    std::getline(std::cin, dummy);
}

void launch_cblas_dgemm(double* A, double* B, double* C){
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, randomMatSize, randomMatSize, randomMatSize, 1.0, A, randomMatSize, B, randomMatSize, 0.1, C, randomMatSize);
}

void FillMatrices(std::vector<std::vector<double>>& matBlock, std::mt19937_64& PRNG, std::uniform_real_distribution<double>& rngdist){
	for(uint32_t i=0; i<3; i++){
		for(uint32_t j=0; j<(randomMatSize*randomMatSize); j++){
			matBlock[i][j] = rngdist(PRNG);
		}
	}
	for(uint32_t i=3; i<(numConcurrentThreads*3); i+=3){
			matBlock[i] = matBlock[0];
			matBlock[i+1] = matBlock[1];
			matBlock[i+2] = matBlock[2];
	}
}

std::mt19937_64 InitPRNG(){
	std::random_device rd;
	std::mt19937_64 PRNG(rd()); //seed PRNG using /dev/urandom or similar OS provided RNG
	std::uniform_real_distribution<double> rngdist{-1.0, 1.0};
	//make sure the internal state of the PRNG is properly mixed by generating 10M random numbers
	//PRNGs often have unreliable distribution uniformity and other statistical properties before their internal state is sufficiently mixed
	for (uint32_t i=0;i<10000000;i++) rngdist(PRNG);
	return PRNG;
}

void PrintMatrices(const std::vector<std::vector<double>>& matBlock){
	for (uint32_t i=0;i<numConcurrentThreads*3;i++){
		std::cout<<i<<std::endl;
		for (uint32_t j=0;j<randomMatSize;j++){
			for (uint32_t k=0;k<randomMatSize;k++){
				std::cout<<matBlock[i][j*randomMatSize + k]<<"  ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
}

int main(){
	std::uniform_real_distribution<double> rngdist{-1.0, 1.0};
	std::vector<std::vector<double>> matBlock(numConcurrentThreads*3);
	std::vector<std::future<void>> futureBlock(numConcurrentThreads);
	
	std::cout<<"*----------------------------*\n";
	std::cout<<"| DGEMM thread safety tester |\n";
	std::cout<<"*----------------------------*\n";
	std::cout<<"Size of random matrices(N=M=K): "<<randomMatSize<<'\n';
	std::cout<<"Number of concurrent calls into OpenBLAS : "<<numConcurrentThreads<<'\n';
	std::cout<<"Number of testing rounds : "<<numTestRounds<<'\n';
	std::cout<<"This test will need "<<(static_cast<uint64_t>(randomMatSize*randomMatSize)*numConcurrentThreads*3*8)/static_cast<double>(1024*1024)<<" MiB of RAM\n"<<std::endl;
	
	std::cout<<"Initializing random number generator..."<<std::flush;
	std::mt19937_64 PRNG = InitPRNG();
	std::cout<<"done\n";
	
	std::cout<<"Preparing to test CBLAS DGEMM thread safety\n";
	std::cout<<"Allocating matrices..."<<std::flush;
	for(uint32_t i=0; i<(numConcurrentThreads*3); i++){
		matBlock[i].resize(randomMatSize*randomMatSize);
	}
	std::cout<<"done\n";
	//pauser();
	std::cout<<"Filling matrices with random numbers..."<<std::flush;
	FillMatrices(matBlock, PRNG, rngdist);
	//PrintMatrices(matBlock);
	std::cout<<"done\n";
	std::cout<<"Testing CBLAS DGEMM thread safety\n";
	omp_set_num_threads(numConcurrentThreads);
	for(uint32_t R=0; R<numTestRounds; R++){
		std::cout<<"DGEMM round #"<<R<<std::endl;
		std::cout<<"Launching "<<numConcurrentThreads<<" threads..."<<std::flush;
		#pragma omp parallel for default(none) shared(futureBlock, matBlock)
		for(uint32_t i=0; i<numConcurrentThreads; i++){
			futureBlock[i] = std::async(std::launch::async, launch_cblas_dgemm, &matBlock[i*3][0], &matBlock[i*3+1][0], &matBlock[i*3+2][0]);
			//launch_cblas_dgemm( &matBlock[i][0], &matBlock[i+1][0], &matBlock[i+2][0]);
		}
		std::cout<<"done\n";
		std::cout<<"Waiting for threads to finish..."<<std::flush;
		for(uint32_t i=0; i<numConcurrentThreads; i++){
			futureBlock[i].get();
		}
		std::cout<<"done\n";
		//PrintMatrices(matBlock);
		std::cout<<"Comparing results from different threads..."<<std::flush;
		for(uint32_t i=3; i<(numConcurrentThreads*3); i+=3){
			for(uint32_t j=0; j<(randomMatSize*randomMatSize); j++){
				if (std::abs(matBlock[i+2][j] - matBlock[2][j]) > 1.0E-13){
					std::cout<<"ERROR: one of the threads returned a different result!"<<i+2<<std::endl;
					std::cout<<"CBLAS DGEMM thread safety test FAILED!"<<std::endl;
					return -1;
				}
			}
		}
		std::cout<<"OK!"<<std::endl;
	}
	std::cout<<"CBLAS DGEMM thread safety test PASSED!"<<std::endl;
	return 0;
}
