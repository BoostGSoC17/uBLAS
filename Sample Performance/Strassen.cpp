#include<bits/stdc++.h>

#define pii pair<int,int>
#define uint unsigned int
#define VT std::vector<T>
#define VVT std::vector<VT>
using namespace std::chrono;

int var = 0;

template<typename T>
void Matrix_Multiply(uint N,T &A, T &B, T &C){
    for(uint i=0;i<N;i++){
        for(uint j=0;j<N;j++){
            C[i][j] = 0;
            for(uint k=0; k<N; k++){
                C[i][j] = C[i][j] + A[i][k]*B[k][j];
            }
        }
    }
}

template<typename T>
void Matrix_Add(int N, VVT &X, VVT &Y, VVT &Z){
    for(uint i=0; i<N; i++){
        for(uint j=0; j<N; j++){
            Z[i][j] = X[i][j] + Y[i][j];
        }
    }
}

template<typename T>
void Matrix_Sub(int N, VVT &X, VVT &Y, VVT &Z){
    for(uint i=0; i<N; i++){
        for(uint j=0; j<N; j++){
            Z[i][j] = X[i][j] - Y[i][j];
        }
    }
}

template<typename T>
void Strassen(uint N, VVT &A, VVT &B, VVT &C){

    if(N == 512) {
        Matrix_Multiply(N, A, B, C);
        return;
    }
    var += 1;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    VVT A11(N,VT(N)), A12(N,VT(N)), A21(N,VT(N)), A22(N,VT(N));
    VVT B11(N,VT(N)), B12(N,VT(N)), B21(N,VT(N)), B22(N,VT(N));
    VVT C11(N,VT(N)), C12(N,VT(N)), C21(N,VT(N)), C22(N,VT(N));
    VVT P1(N,VT(N)), P2(N,VT(N)), P3(N,VT(N)), P4(N,VT(N)),P5(N,VT(N)), P6(N,VT(N)), P7(N,VT(N));
    VVT AA(N,VT(N)), BB(N,VT(N));

    uint mid = N>>1;
    for(uint i=0; i<N; i++){
        for(uint j=0;j<N;j++){
            if(i<mid && j<mid)
                A11[i][j] = A[i][j];
            else if(i<mid)
                A12[i][j-mid] = A[i][j];
            else if(i>=mid && j<mid)
                A21[i-mid][j] = A[i][j];
            else
                A22[i-mid][j-mid] = A[i][j];

            if(i<mid && j<mid)
                B11[i][j] = B[i][j];
            else if(i<mid)
                B12[i][j-mid] = B[i][j];
            else if(i>=mid && j<mid)
                B21[i-mid][j] = B[i][j];
            else
                B22[i-mid][j-mid] = B[i][j];
        }
    }

    Matrix_Add(N>>1, A11, A22, AA);
    Matrix_Add(N>>1, B11, B22, BB);
    Strassen(N>>1, AA, BB, P1);

    Matrix_Add(N>>1, A21, A22, AA);
    Strassen(N>>1, AA, B11, P2);

    Matrix_Sub(N>>1, B12, B22, BB);
    Strassen(N>>1, A11, BB, P3);

    //Calculate M4 = A3 × (B2 - B0)
    Matrix_Sub((N>>1), B21, B11, BB);
    Strassen((N>>1), A22, BB, P4);

    //Calculate M5 = (A0 + A1) × B3
    Matrix_Add((N>>1), A11, A12, AA);
    Strassen((N>>1), AA, B22, P5);

    //Calculate M6 = (A2 - A0) × (B0 + B1)
    Matrix_Sub((N>>1), A21, A11, AA);
    Matrix_Add((N>>1), B11, B12, BB);
    Strassen((N>>1), AA, BB, P6);

    //Calculate M7 = (A1 - A3) × (B2 + B3)
    Matrix_Sub((N>>1), A12, A22, AA);
    Matrix_Add((N>>1), B21, B22, BB);
    Strassen((N>>1), AA, BB, P7);

    //Calculate C0 = M1 + M4 - M5 + M7
    Matrix_Add((N>>1), P1, P4, AA);
    Matrix_Sub((N>>1), P7, P5, BB);
    Matrix_Add((N>>1), AA, BB, C11);

    //Calculate C1 = M3 + M5
    Matrix_Add((N>>1), P3, P5, C12);

    //Calculate C2 = M2 + M4
    Matrix_Add((N>>1), P2, P4, C21);

    //Calculate C3 = M1 - M2 + M3 + M6
    Matrix_Sub((N>>1), P1, P2, AA);
    Matrix_Add((N>>1), P3, P6, BB);
    Matrix_Add((N>>1), AA, BB, C22);

    //Set the result to C[][N]
    for(int i=0; i<N; i++) {
       for(int j=0; j<N; j++) {
          if(i<mid && j<mid)
            C[i][j] = C11[i][j];
          else if(i<mid)
            C[i][j] = C12[i][j-mid];
          else if(i>=mid && j<mid)
            C[i][j] = C21[i-mid][j];
          else
            C[i][j] = C22[i-mid][j-mid];
       }
    }

    A11.clear(); A12.clear(); A21.clear(); A22.clear();
    B11.clear(); B12.clear(); B21.clear(); B22.clear();
    C11.clear(); C12.clear(); C21.clear(); C22.clear();
    P1.clear(); P2.clear(); P3.clear(); P4.clear(); P5.clear(); P6.clear(); P7.clear();
    AA.clear(); BB.clear();

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << std::fixed;
    std::cout << std::setprecision(7) << time_span.count() << " " << var << "\n";

    var -= 1;
}

template<typename T>
void multiply(uint N,VVT &A, VVT &B, VVT &C){
    if(N<512){
        Matrix_Multiply(N,A,B,C);
        return;
    }
    Strassen(N,A,B,C);
}

const uint N = 1024;

int main(){

    std::vector<std::vector<int> > A(N,std::vector<int>(N,1));
    std::vector<std::vector<int> > B(N,std::vector<int>(N,1));
    std::vector<std::vector<int> > C(N,std::vector<int>(N,1));

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    multiply(N,A,B,C);
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
