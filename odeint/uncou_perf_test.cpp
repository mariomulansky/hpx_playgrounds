#include <iostream>
#include <vector>
#include <cmath>

#include <boost/numeric/odeint.hpp>

typedef std::vector< double > dvec;

void( const dvec &x , dvec &dxdt , double )
{
    for( int i=0 ; i<N ; ++i )
    {
        dxdt[i] = pow(abs(x[i]),0.5) * sin(x[i]) + cos(x[i]);
    }
}

const size_t M = 10000000;
const size_t N = 4 * M;

int main()
{
    dvec x( N , 100.0 );
    
}
