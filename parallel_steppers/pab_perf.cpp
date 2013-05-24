// Copyright Mario Mulansky 2013

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <cmath>

#define HPX_LIMIT 10

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/numeric/odeint/stepper/parallel_extrapolation_stepper.hpp>
#include <boost/numeric/odeint/stepper/parallel_adams_bashforth_stepper.hpp>

#include <boost/foreach.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp> 

#include <boost/math/special_functions/sign.hpp>

#include "dataflow_shared_wrapper.hpp"
#include "dataflow_shared_algebra.hpp"
#include "dataflow_operations.hpp"
#include "dataflow_shared_copy.hpp"
#include "identity_action.hpp"
#include "serialization.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

using boost::numeric::odeint::parallel_extrapolation_stepper;
using boost::numeric::odeint::parallel_adams_bashforth_stepper;

const size_t N = 10000;
const double GAMMA = 1.2;

typedef std::array< double , N > darray;
typedef std::shared_ptr< darray > shared_array;
typedef dataflow_base< shared_array > state_type;

typedef parallel_extrapolation_stepper< 8 , 
                                        state_type , double , state_type , double , 
                                        dataflow_shared_algebra , 
                                        dataflow_operations 
                                        > extrapolation_stepper_type;

typedef parallel_adams_bashforth_stepper< 3 , 
                                        state_type , double , state_type , double , 
                                        dataflow_shared_algebra , 
                                        dataflow_operations 
                                        > pab_stepper_type;


inline double coupling( const double x )
{
    using std::pow;
    using std::abs;
    using boost::math::sign;
    //return pow( abs(sin( x )) , 1.35 ) * sign(sin(x)) - GAMMA * ( 1.0 - cos( x ) );
    return sin( x ) - GAMMA * ( 1.0 - cos( x ) );
}

shared_array rhs_func( const shared_array x_ , shared_array dxdt_ )
{
    // hpx::cout << "rhs start\n" << hpx::flush;
    darray& x = *x_;
    darray& dxdt = *dxdt_;
    dxdt[0] = coupling( x[1]-x[0] );
    for( size_t i=1 ; i<N-1 ; ++i )
    {
        dxdt[i] = coupling( x[i+1]-x[i] ) + coupling( x[i-1]-x[i] );
    }
    dxdt[N-1] = coupling( x[N-2] - x[N-1] );
    // hpx::cout << "rhs end\n" << hpx::flush;
    return dxdt_;
}
HPX_PLAIN_ACTION( rhs_func , rhs_func_action )

void rhs( const state_type &x , state_type &dxdt , double t )
{
    dxdt = dataflow< rhs_func_action >( find_here() , x , dxdt );
}

void sync_swap( state_type &x1 , state_type &x2 , state_type &sync1 , state_type &sync2 )
{
    state_type tmp = x1;
    typedef identity_sync3_action< shared_array , shared_array > sync_ident3;
    x1 = dataflow< sync_ident3 >( find_here() , x2 , x1 , sync1 , sync2 );
    typedef identity_sync1_action< shared_array , shared_array > sync_ident1;
    x2 = dataflow< sync_ident1 >( find_here() , tmp , x1 );
    sync1 = dataflow< sync_ident1 >( find_here() , sync1 , x2 );
    sync2 = dataflow< sync_ident1 >( find_here() , sync2 , x2 );
}

int hpx_main(boost::program_options::variables_map& vm)
{
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const std::size_t M = vm["M"].as<std::size_t>();
    const double dt = 0.1;

    double mean_time = 0.0;

    for( size_t m=0 ; m<M ; ++m )
    {
        darray in;
        std::uniform_real_distribution<double> distribution( 0.0 , 2*3.14159 );
        std::mt19937 engine( 0 ); // Mersenne twister MT19937
        auto generator = std::bind(distribution, engine);
        std::generate( in.begin() , in.end() , std::ref(generator) );
        typedef initialize_shared_action<darray> init_action;
        state_type x_in = dataflow< init_action >( find_here() , in );
        darray out = {{0.0}};
        state_type x_out = dataflow< init_action >( find_here() , out );

        pab_stepper_type stepper;
        // for some reason this is necessary
        wait( x_in.get_future() );
        wait( x_out.get_future() );
        //hpx::cout << (boost::format("%f\n") % (*(x_in.get_future().get()))[0]) << hpx::flush;
        hpx::util::high_resolution_timer timer;
        
        for( size_t n = 0 ; n < steps ; ++n )
        {
            stepper.do_step( rhs , x_in , n*dt , x_out , 0.1 );
            // swap in <-> out and synchronize internal states of the stepper
            sync_swap( x_in , x_out , stepper.m_states[0].m_v , stepper.m_states[1].m_v );
        }
        wait( x_in.get_future() );

        mean_time += timer.elapsed();

        std::clog << (boost::format("runtime: %fs\n") %timer.elapsed());

        std::clog << (boost::format("%f\n") % (*(x_in.get_future().get()))[0]);
    }
    
    hpx::cout << (boost::format("mean runtime: %f\n") % (mean_time/M) ) << hpx::flush;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "steps",
          boost::program_options::value<std::size_t>()->default_value(100),
          "Steps (100)")
        ;

    desc_commandline.add_options()
        ( "M",
          boost::program_options::value<std::size_t>()->default_value(10),
          "Runs (10)")
        ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
