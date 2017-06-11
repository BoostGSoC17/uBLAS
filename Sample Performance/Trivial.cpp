#include<bits/stdc++.h>

#define uint unsigned int
#define VT std::vector<T>
#define VVT std::vector<VT>
using namespace std::chrono;

const uint N = 2048;

template<typename T>
void Multiply(uint N, VVT &A, VVT &B, VVT &C){
    for(uint i=0;i<N;i++){
        for(uint j=0;j<N;j++){
            C[i][j] = 0;
            for(uint k=0;k<N;k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main(){

    std::vector<std::vector<int> > A(N,std::vector<int>(N,1));
    std::vector<std::vector<int> > B(N,std::vector<int>(N,1));
    std::vector<std::vector<int> > C(N,std::vector<int>(N,1));

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    Multiply(N,A,B,C);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << std::fixed;
    std::cout << std::setprecision(7) << time_span.count() << "\n";

    /*for(uint i=0;i<C.size();i++){
        for(uint j=0;j<C.size();j++){
            cout << C[i][j] << "\t";
        }
        cout << "\n";
    }*/

    return 0;
}
