/* Stress test for library shutdown racing with in-flight BLAS calls
 * (https://github.com/OpenMathLib/OpenBLAS/issues/5954).
 *
 * Windows only. On POSIX, exit() runs the library destructor while worker
 * threads are still computing into OpenBLAS-owned buffers, which no amount of
 * locking inside blas_shutdown can make safe, so there is nothing to assert
 * there; CMakeLists.txt only registers this test on WIN32.
 *
 * The parent re-executes itself as short-lived children and checks that each
 * one terminates cleanly, turning shutdown-path crashes and deadlocks into
 * ordinary test failures. Each child (--child-storm N) starts N callers that
 * allocate their matrices and park on a gate, releases them so they all enter
 * their first dgemm at once, and exits a millisecond later while that
 * allocation storm is still in flight.
 *
 * N must exceed NUM_BUFFERS = MAX(50, NUM_THREADS * 2 * NUM_PARALLEL) for the
 * build under test; below that every slot is already mapped and the race is
 * unreachable.
 */
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#ifdef OPENBLAS_USE_GENERATED_CBLAS_H
#include "generated/cblas.h"
#else
#include "../cblas.h"
#endif

#include <windows.h>

namespace {

const blasint stormM = 200, stormK = 120, stormN = 90; /* the gh-5954 shape */
const blasint poolDim = 320; /* above the multithreading threshold, so the pool spins up */
const uint32_t defaultStormCallers = 128;
const uint32_t stormDelayMs = 3; /* gate to sweep; at 0 the sweep beats the allocations */
const int stormBlasThreads = 4;
const int stormTimeoutSec = 15;
const int numStormChildren = 40;

std::atomic<uint32_t> parked(0); /* callers built and waiting on the gate */
std::atomic<bool> gate(false);

void fillOperands(std::vector<double>& A, std::vector<double>& B) {
	for (size_t i = 0; i < A.size(); i++) A[i] = (i % 1000) / 1000.0;
	for (size_t i = 0; i < B.size(); i++) B[i] = (i % 997) / 997.0;
}

void dgemmOnce(blasint m, blasint k, blasint n) {
	std::vector<double> A(m * k), B(k * n), C(m * n);
	fillOperands(A, B);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
	            1.0, A.data(), m, B.data(), k, 0.1, C.data(), m);
}

/* Allocate before parking, so that when the gate opens nothing stands between
   the thread and its first dgemm. */
void gatedWorker(blasint m, blasint k, blasint n) {
	std::vector<double> A(m * k), B(k * n), C(m * n);
	fillOperands(A, B);
	parked.fetch_add(1, std::memory_order_release);
	while (!gate.load(std::memory_order_acquire)) std::this_thread::yield();
	for (;;)
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
		            1.0, A.data(), m, B.data(), k, 0.1, C.data(), m);
}

int ChildStorm(uint32_t nCallers) {
	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
	openblas_set_num_threads(stormBlasThreads);

	/* Build the OpenBLAS worker pool first, so the storm is buffer allocation
	   and not pool startup. */
	dgemmOnce(poolDim, poolDim, poolDim);

	for (uint32_t i = 0; i < nCallers; i++)
		std::thread(gatedWorker, stormM, stormK, stormN).detach();
	for (int ms = 0; parked.load(std::memory_order_acquire) < nCallers && ms < 10000; ms++)
		std::this_thread::sleep_for(std::chrono::milliseconds(1));

	gate.store(true, std::memory_order_release);
	std::this_thread::sleep_for(std::chrono::milliseconds(stormDelayMs));
	std::exit(0);
}

/* Returns 0 if the child exited cleanly, nonzero otherwise; fills outcome. */
int RunChild(const std::string& args, int timeoutSec, std::string& outcome) {
	char exe[MAX_PATH];
	if (GetModuleFileNameA(NULL, exe, MAX_PATH) == 0) {
		outcome = "GetModuleFileName failed";
		return 1;
	}
	std::string cmd = "\"" + std::string(exe) + "\" " + args;

	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	if (!CreateProcessA(NULL, &cmd[0], NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
		outcome = "CreateProcess failed";
		return 1;
	}

	int ret = 1;
	char buf[64];
	if (WaitForSingleObject(pi.hProcess, timeoutSec * 1000) != WAIT_OBJECT_0) {
		TerminateProcess(pi.hProcess, 1);
		WaitForSingleObject(pi.hProcess, 5000);
		snprintf(buf, sizeof(buf), "HANG (killed after %ds)", timeoutSec);
	} else {
		DWORD code = 1;
		GetExitCodeProcess(pi.hProcess, &code);
		if (code == 0) {
			snprintf(buf, sizeof(buf), "clean exit");
			ret = 0;
		} else {
			snprintf(buf, sizeof(buf), "CRASH (exit code 0x%08lX)", (unsigned long)code);
		}
	}
	outcome = buf;

	CloseHandle(pi.hThread);
	CloseHandle(pi.hProcess);
	return ret;
}

} // namespace

int main(int argc, char* argv[]) {
	if (argc >= 3 && std::strcmp(argv[1], "--child-storm") == 0)
		return ChildStorm(uint32_t(std::atoi(argv[2])));
	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);

	uint32_t callers = defaultStormCallers;
	if (argc >= 2) {
		int n = std::atoi(argv[1]);
		if (n > 0) callers = uint32_t(n);
	}

	int failures = 0;
	std::cout << "Testing process exit during an allocation storm (" << callers << " callers)"
	          << std::endl;
	for (int i = 0; i < numStormChildren; i++) {
		std::string outcome;
		failures += RunChild("--child-storm " + std::to_string(callers), stormTimeoutSec, outcome);
		std::cout << "  storm child " << i << ": " << outcome << std::endl;
	}

	if (failures) {
		std::cout << "CBLAS DGEMM shutdown safety test FAILED! (" << failures
		          << " child processes)" << std::endl;
		return 1;
	}
	std::cout << "CBLAS DGEMM shutdown safety test PASSED!" << std::endl;
	return 0;
}
