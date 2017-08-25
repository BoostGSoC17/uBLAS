#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/timer.hpp>

using namespace boost::numeric::ublas;

int main() {
	boost::timer timer;
	double elapsed_time;

	matrix<int> A, B, C, D, E;
	int AR = 100, BR = 200, CR = 300, DR = 400, DC = 500;
	while(DC <= 2000) {
		A.resize(AR,BR); B.resize(BR,CR); C.resize(CR,DR); D.resize(DR,DC);

		for(int i=0;i<A.size1();i++) {
			for(int j=0;j<A.size2();j++){
				A(i,j) = 1;
			}
		}
		for(int i=0;i<B.size1();i++) {
			for(int j=0;j<B.size2();j++){
				B(i,j) = 1;
			}
		}
		for(int i=0;i<C.size1();i++) {
			for(int j=0;j<C.size2();j++){
				C(i,j) = 1;
			}
		}
		for(int i=0;i<D.size1();i++) {
			for(int j=0;j<D.size2();j++){
				D(i,j) = 1;
			}
		}

		timer.restart();
		E = A * B * C * D;
		elapsed_time = timer.elapsed();
		std::cout << "(" << AR << ", " << BR << ") * (";
		std::cout << "(" << BR << ", " << CR << ") * (";
		std::cout << "(" << CR << ", " << DR << ") * (";
		std::cout << "(" << DR << ", " << DC << ")\n";
		std::cout << "Time taken:\t" << elapsed_time << "\n";
		AR += 100; BR += 100; CR += 100; DR += 100; DC += 100;
	}

	return 0;
}