#include <bits/stdc++.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

using namespace std::chrono;
namespace ublas = boost::numeric::ublas;

int main () {
	ublas::matrix<int> A, B, ans;
	//freopen("Benchmark.txt", "w", stdout);
    int N = 100;
	while(N<=2000) {
		A.resize(N,N+50); B.resize(N+50,N), ans.resize(N,N);
		for(int i=0;i<A.size1();i++) {
			for(int j=0;j<A.size2();j++) {
				A(i,j) = 1; B(i,j) = 1;			
			}		
		}
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		ans = prod(A, B);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		
		//std::cout << "Size = " << N << "\t Time = ";
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    	std::cout << std::fixed;
    	std::cout << std::setprecision(7) << time_span.count() << "\n";	
		N += 100;
	}
	
	
    return 0;
}
