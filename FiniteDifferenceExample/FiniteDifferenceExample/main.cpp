#include "FiniteDifferenceExample.h"

void test(int N, std::string mesh, double r_e, std::string convection, std::string path) {
	FiniteDifferenceExample FDE(1.0, 1.0, 0.02, 1.0, 0.0, 1.0, N, mesh, r_e, convection);
	Data data;
	FDE.solve(data);
	data.save(path, convection);
}

void test1() {
	test(10, "uniform", 0.0, "CDS", "test1");
	test(10, "uniform", 0.0, "UDS", "test1");
}

void test2() {
	test(40, "uniform", 0.0, "CDS", "test2");
	test(40, "uniform", 0.0, "UDS", "test2");
}

void test3() {
	test(10, "non-uniform", 0.7, "CDS", "test3");
	test(10, "non-uniform", 0.7, "UDS", "test3");
}

int main() {
	test1();
	test2();
	test3();
	return 0;
}
