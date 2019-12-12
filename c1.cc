
#include <iostream>
#include <cmath>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/make-candidate.hpp>
#include <kv/psa.hpp>
#include <kv/autodif.hpp>
#include <kv/ode-param.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;
template <class T, class F>
// return a subset int b
bool enclosure(ub::vector<interval<T> > a, ub::vector<interval<T> > b){
    if(a.size() != b.size())
        return false;
    else{
        int n = a.size(); 
        while(n-- > 0){
            if(a[n].lower() < b[n].lower() || b[n].upper() < a[n].upper())
                return false;
        }
    }
    return true;
}
template <class T, class F>
T LogarithmicNorm(ub::matrix<T> &q){
    if(q.size1() != q.size2())
        return -1;
    T mu = 0;
    for(int i = 0; i < q.size1(); i++){
        T tmp = q[i][i];
        for(int k = 0; k < q.size1(); k++){
            if(k == i)
                continue;
            tmp += q[k][i];
        }
        if(mu < tmp)
            mu = tmp;

    }

    return mu;
}
template <class T, class F>
int 
ode_c1lohner(F f, F dfdx, ub::vector< interval<T> >& init_x,ub::matrix< interval<T> >& init_v, const T& start, T& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL) {
int n = init_x.size();

ub::vector<interval<T> > W1(n), W2(n), Y(n), tmp(n);
ub::matrix<interval<T> > W3(n,n);
ub::identity_matrix<T> I(n, n);
T h = end - start;
//part1
Y = tmp = I;
T max_ratio = std::max(max_ratio, width(Y) / width(Y));
while (!enclosure(tmp, Y))
{
    Y = tmp;
    tmp = init_x + interval<T>(0, h) * f(Y, start);
}
W1 = intersect(Y, tmp);

//part2
T l = LogarithmicNorm( dfdx( W1, start ) );
ub::matrix<interval<T> > W(n, n);
for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
        W[i][j] = interval<T>(-exp(l * h), exp(l * h));
W3 = intersect(W, I + interval<T>(0, h) * dfdx(W1, start) * W);
//part3
ub::matrix<interval<T>> A(n, n), J(n, n);
//part4

return 0;
}
}
int main(void){
    return 0;
}
