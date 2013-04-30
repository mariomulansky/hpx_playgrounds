// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>

#include <sys/time.h>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>

#include "../odeint_shared/dataflow_shared_resize.hpp"
#include "../odeint_shared/dataflow_shared_algebra.hpp"
#include "../odeint_shared/dataflow_shared_operations.hpp"
#include "../odeint_shared/serialization.hpp"
#include "mixed_toda_nl_system.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

using hpx::cout;
using hpx::flush;

using boost::format;
using boost::io::group;

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;

typedef dataflow_base< double > df_base;
typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef std::vector< dataflow_base< shared_vec > > state_type;
typedef std::vector< dvec > vec_dvec;

typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       dataflow_shared_algebra ,
                                       dataflow_operations > stepper_type;


shared_vec sync_identity3( shared_vec x , shared_vec sync1 , shared_vec sync2 , shared_vec sync3 )
{
    return x;
}

HPX_PLAIN_ACTION( sync_identity3 , sync_identity3_action )

shared_vec sync_identity2( shared_vec x , shared_vec sync1 , shared_vec sync2 )
{
    return x;
}

HPX_PLAIN_ACTION( sync_identity2 , sync_identity2_action )

shared_vec sync_identity1( shared_vec x , shared_vec sync1 )
{
    return x;
}

HPX_PLAIN_ACTION( sync_identity1 , sync_identity1_action )

// syncronized swap emulating nearest neighbor coupling
void synchronized_swap( state_type &x_in , state_type &x_out )
{
    dataflow_base< shared_vec > x_left = x_out[0];
    dataflow_base< shared_vec > x_tmp = dataflow< sync_identity2_action >( find_here() , x_in[0] , 
                                                                           x_out[0] , x_out[1] );
    x_in[0] = dataflow< sync_identity1_action >( find_here() , x_out[0] , x_tmp );
    x_out[0] = x_tmp;
    x_tmp = dataflow< sync_identity3_action >( find_here() , x_in[1] , 
                                               x_out[1] , x_left , x_out[2] );
    x_left = x_out[1];
    x_in[1] = dataflow< sync_identity1_action >( find_here() , x_out[1] , x_tmp );
    x_out[1] = x_tmp;
    x_tmp = dataflow< sync_identity2_action >( find_here() , x_in[2] , 
                                               x_out[2] , x_left );
    x_in[2] = dataflow< sync_identity1_action >( find_here() , x_out[2] , x_tmp );
    x_out[2] = x_tmp;
}


int hpx_main(boost::program_options::variables_map& vm)
{
    const std::size_t M = vm["M"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();

    const std::size_t N = 3;

    hpx::cout << boost::format( "number of dataflows: %d, number of elements per dataflow: %d, steps: %d, dt: %.3f\n") % N % M % steps % dt ;

    state_type q_in( N );
    state_type p_in( N );
    state_type q_out( N );
    state_type p_out( N );

    for( size_t i=0 ; i<N ; ++i )
    {
        q_in[i] = dataflow< initialize_action >( find_here() , 
                                                 std::allocate_shared< dvec >( std::allocator<dvec>() ) ,
                                                 M ,
                                                 0.0 );
        p_in[i] = dataflow< initialize_random_action >( find_here() , 
                                                        std::allocate_shared< dvec >( std::allocator<dvec>() ) ,
                                                        M ,
                                                        i );
        q_out[i] = dataflow< initialize_action >( find_here() , 
                                                  std::allocate_shared< dvec >( std::allocator<dvec>() ) ,
                                                  M ,
                                                  0.0 );
        p_out[i] = dataflow< initialize_action >( find_here() , 
                                                  std::allocate_shared< dvec >( std::allocator<dvec>() ) ,
                                                  M ,
                                                  0.0 );

    }

    std::vector< future<shared_vec> > futures_q( N );
    std::vector< future<shared_vec> > futures_p( N );
    for( size_t i=0 ; i<N ; ++i )
    {
        futures_q[i] = q_in[i].get_future();
        futures_p[i] = p_in[i].get_future();
    }

    wait( futures_q );
    wait( futures_p );
    hpx::cout << boost::format("Initialization complete, energy: %.2f\n") % (energy( futures_q , futures_p ));

    hpx::util::high_resolution_timer timer;
    
    // do the numerical integration using dataflow objects
    stepper_type stepper;
    for( size_t t=0 ; t<steps ; ++t )
    {
        auto in = std::make_pair( boost::ref(q_in) , boost::ref(p_in) );
        auto out = std::make_pair( boost::ref(q_out) , boost::ref(p_out) );
        stepper.do_step( mixed_toda_nl_system , 
                         in ,
                         t*dt , 
                         out , 
                         dt );
        synchronized_swap( q_in , q_out );
        synchronized_swap( p_in , p_out );
    }

    for( size_t i=0 ; i<N ; ++i )
    {
        futures_q[i] = q_in[i].get_future();
        futures_p[i] = p_in[i].get_future();
    }

    // here we wait for the results
    wait( futures_q );
    wait( futures_p );
    //wait( futures_out );

    cout << (boost::format("%d") % (M) );
    cout << '\t' << (boost::format("%f") % (timer.elapsed())) << "\n";
    
    cout << boost::format( "Final energy: %.2f\n" ) % energy( futures_q , futures_p ) << flush;

    // print the values
    // for( size_t i=0 ; i<N ; ++i )
    //     for( size_t m=0 ; m<M ; m++ )
    //     {
    //         std::cout << (*(futures_q[i].get()))[m] << ',';
    //         std::cout << (*(futures_p[i].get()))[m] << '\t';
    //     }
    // std::cout << std::endl;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "M",
          boost::program_options::value<std::size_t>()->default_value(128),
          "M - Number of elements in each part (128)")
        ;
    desc_commandline.add_options()
        ( "steps",
          boost::program_options::value<std::size_t>()->default_value(100),
          "time steps (100)")
        ;
    desc_commandline.add_options()
        ( "dt",
          boost::program_options::value<double>()->default_value(0.01),
          "step size (0.01)")
        ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
