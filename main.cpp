#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

int main () {
    using namespace boost::numeric::ublas;
    matrix<int> A (5, 3), X(3,3);
    matrix<int> B(3,4), C(4,3), tmp(3,3);
    for (unsigned i = 0; i < A.size1 (); ++ i)
        for (unsigned j = 0; j < A.size2 (); ++ j)
            A(i,j) = 1;

    for (unsigned i = 0; i < X.size1 (); ++ i)
        for (unsigned j = 0; j < X.size2 (); ++ j)
            X(i,j) = 1;

    for (unsigned i = 0; i < B.size1 (); ++ i)
        for (unsigned j = 0; j < B.size2 (); ++ j)
            B(i,j) = 1;

    for (unsigned i = 0; i < C.size1 (); ++ i)
        for (unsigned j = 0; j < C.size2 (); ++ j)
            C(i,j) = 1;

    //tmp = prod(A, X - matrix<int>(prod(B,C)));
    //tmp = prod(matrix<int>(prod(B,C)), X);
    tmp = prod(A,B);
//    tmp = A + X;
    std::cout << A << std::endl << X << std::endl << tmp << std::endl;// << std::endl << tmp;
}
