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

using boost::timer::cpu_timer;
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
    int N2 = 1024;
    int init_length = 128;
    int steps = 10;
    double dt = 0.1;
    if( argc > 1 )
        N1 = atoi( argv[1] );
    if( argc > 2 )
        N2 = atoi( argv[2] );
    int block_size = 8;
    if( argc > 3 )
        steps = atoi( argv[3] );

    std::clog << "Size: " << N1 << "x" << N2 << " with " << steps << " steps" << std::endl;

    while( block_size <= N1/8 )
    {

        double avrg_time = 0;

        for( size_t n=0 ; n<10 ; ++n )
        {

            lattice2d system( KAPPA , LAMBDA , beta , block_size );

            // initialize
            state_type q( N1 , dvec( N2 , 0.0 ) );
            state_type p( N1 , dvec( N2 , 0.0 ) );
    
            //fully random
            for( size_t i=0 ; i<N1 ; ++i )
            {
                std::uniform_real_distribution<double> distribution( 0.0 );
                std::mt19937 engine( i ); // Mersenne twister MT19937
                auto generator = std::bind( distribution , engine );
                std::generate( p[i].begin() , p[i].end() , generator );
            }

            // std::clog << "# Initial energy: " << system.energy( q , p ) << " (fully random)" << std::endl;

            //partly random
            // for( size_t i=N1/2-init_length/2 ; i<N1/2+init_length/2 ; ++i )
            // {
            //     std::uniform_real_distribution<double> distribution( 0.0 );
            //     std::mt19937 engine( i ); // Mersenne twister MT19937
            //     auto generator = std::bind( distribution , engine );
            //     std::generate( p[i].begin()+N2/2-init_length/2 ,
            //                    p[i].begin()+N2/2+init_length/2 ,
            //                    generator );
            // }

            // for( int i=0 ; i<N1 ; ++i )
            // {
            //     for( int j=0 ; j<N2 ; ++j )
            //     {
            //         std::cout << q[i][j] << "," << p[i][j] << '\t';
            //     }
            //     std::cout << std::endl;
            // }

            spreading_observer obs( KAPPA , LAMBDA , beta , block_size );

            //std::cout << "# Initial energy: " << system.energy( q , p ) << std::endl;
    
            cpu_timer timer;

            integrate_n_steps( stepper_type( nested_omp_algebra<range_algebra>( block_size ) ) , 
                               system , 
                               std::make_pair( std::ref(q) , std::ref(p) ) , 
                               0.0 , dt , steps );
            avrg_time += static_cast<double>(timer.elapsed().wall)/(1000*1000*1000);
            // std::ref(obs) );
            // std::clog << "Final energy: " << system.energy( q , p ) << std::endl;
        }
    
        std::cout << block_size << '\t' << avrg_time/(10) << std::endl;

        block_size *= 2;
    }

    return 0;
}
