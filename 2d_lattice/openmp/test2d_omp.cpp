// Copyright 2013 Mario Mulansky

#include <iostream>
#include <vector>
#include <random>

#include <boost/numeric/odeint.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include "lattice2d.hpp"
#include "nested_range_algebra_omp.hpp"
#include "resize.hpp"
#include "spreading_observer.hpp"

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;
using boost::numeric::odeint::range_algebra;

using boost::timer::auto_cpu_timer;
using boost::timer::cpu_times;

typedef std::vector< double > dvec;
typedef std::vector< dvec > state_type;

typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       nested_omp_algebra<range_algebra> > stepper_type;

const double KAPPA = 3.5;
const double LAMBDA = 4.5;
const double beta = 1.0;

int main( int argc , char* argv[] )
{
    int N1 = 1024;
    int N2 = 1024;\
    int init_length = 128;
    int steps = 100;
    double dt = 0.01;
    if( argc > 1 )
        N1 = atoi( argv[1] );
    if( argc > 2 )
        N2 = atoi( argv[2] );
    int block_size = N1/4;
    if( argc > 3 )
        block_size = atoi( argv[3] );
    if( argc > 4 )
        init_length = atoi( argv[4] );    
    if( argc > 5 )
        steps = atoi( argv[5] );
    if( argc > 6 )
        dt = atof( argv[6] );

    std::cout << "Size: " << N1 << "x" << N2 << " with " << block_size << " rows per thread" << std::endl;

    // initialize
    state_type q( N1 , dvec( N2 , 0.0 ) );
    state_type p( N1 , dvec( N2 , 0.0 ) );

    //fully random
    // for( size_t i=0 ; i<N1 ; ++i )
    // {
    //     std::uniform_real_distribution<double> distribution( 0.0 );
    //     std::mt19937 engine( i ); // Mersenne twister MT19937
    //     auto generator = std::bind( distribution , engine );
    //     std::generate( p[i].begin() , p[i].end() , generator );
    // }

    //partly random
    for( size_t i=N1/2-init_length/2 ; i<N1/2+init_length/2 ; ++i )
    {
        std::uniform_real_distribution<double> distribution( 0.0 );
        std::mt19937 engine( i ); // Mersenne twister MT19937
        auto generator = std::bind( distribution , engine );
        std::generate( p[i].begin()+N2/2-init_length/2 ,
                       p[i].begin()+N2/2+init_length/2 ,
                       generator );
    }

    // for( int i=0 ; i<N1 ; ++i )
    // {
    //     for( int j=0 ; j<N2 ; ++j )
    //     {
    //         std::cout << q[i][j] << "," << p[i][j] << '\t';
    //     }
    //     std::cout << std::endl;
    // }


    lattice2d system( KAPPA , LAMBDA , beta , block_size );
    spreading_observer obs( KAPPA , LAMBDA , beta , block_size );

    std::cout << "Initial energy: " << system.energy( q , p ) << std::endl;

    {
    auto_cpu_timer timer( 3 , "%w sec\n");

    integrate_n_steps( stepper_type( nested_omp_algebra<range_algebra>( block_size ) ) , 
                       system , 
                       std::make_pair( std::ref(q) , std::ref(p) ) , 
                       0.0 , dt , steps );
                       //std::ref(obs) );

    }

    // std::cout << "Time: " << elapsed << std::endl;
    std::cout << "Final energy: " << system.energy( q , p ) << std::endl;

    std::pair< double , double > x;
    BOOST_FOREACH( x , obs.m_values )
    {
        std::cout << x.first << "\t" << x.second << std::endl;
    }

    // for( int i=0 ; i<N1 ; ++i )
    // {
    //     for( int j=0 ; j<N2 ; ++j )
    //     {
    //         std::cout << q[i][j] << "," << p[i][j] << '\t';
    //     }
    //     std::cout << std::endl;
    // }
    return 0;
}
