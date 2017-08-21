//
//  Copyright (c) 2017 Mohd Sharique
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/timer.hpp>

using namespace boost::numeric::ublas;

int main() {
	boost::timer timer;
	double elapsed_time;

	int N=200;
	matrix<int> A, B, C;
	while(N <= 2000) {
		A.resize(N,N-100); B.resize(N-100,N);
		for(unsigned int i=0; i<N; i++) {
			for(unsigned int j=0; j<N; j++) {
				A(i,j) = 1;	B(i, j) = 1;			
			} 
		}
		timer.restart();
		C = prod(A, B);
		elapsed_time = timer.elapsed();
		std::cout << "Matrix dimensions:\t(" << N << ", " << N-100 << ")\t(" << N-100 << ", " << N << "\n"; 
		std::cout << "Time taken:\t" << elapsed_time << "\n";
		N += 100;
	}

	return 0;
}

