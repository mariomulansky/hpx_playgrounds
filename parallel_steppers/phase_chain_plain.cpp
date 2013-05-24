#include <iostream>
#include <array>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>

#include <boost/numeric/odeint/stepper/parallel_extrapolation_stepper.hpp>

using boost::numeric::odeint::parallel_extrapolation_stepper;

const size_t N=4096;
const double GAMMA=1.2;

typedef std::array< double , N > state_type;

typedef parallel_extrapolation_stepper< 8 , state_type > stepper_type;

inline double coupling( const double x )
{
    return sin( x ) - GAMMA * ( 1.0 - cos( x ) );
}

void rhs( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = coupling( x[1]-x[0] );
    for( size_t i=1 ; i<N-1 ; ++i )
    {
        dxdt[i] = coupling( x[i+1]-x[i] ) + coupling( x[i-1]-x[i] );
    }
    dxdt[N-1] = coupling( x[N-2] - x[N-1] );
}

int main( int argc , char **argv )
{
    int steps = 100;
    if( argc>1 )
        steps = atoi( argv[1] );
    double dt = 0.1;

    std::cout.precision(12);

    state_type x , x_out;
    std::uniform_real_distribution<double> distribution( 0.0 , 2*3.14159 );
    std::mt19937 engine( 0 ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);
    std::generate( x.begin() , x.end() , std::ref(generator) );

    std::cout << x[0] << std::endl;

    stepper_type stepper;
    double tic = clock(); //TIC
    for( int n=0 ; n<steps ; ++n )
    {
        stepper.do_step( rhs , x , n*dt , dt );
        //x = x_out;
    }
    double toc = clock();   //TOC
    
    std::cout << "It took "<< (toc-tic)/CLOCKS_PER_SEC <<" second(s)."<< std::endl;

    std::cout << x[0] << std::endl;
}
