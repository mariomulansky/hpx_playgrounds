// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/util/unwrap.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/foreach.hpp>

#include "local_dataflow_shared_resize.hpp"
#include "local_dataflow_algebra.hpp"
#include "local_dataflow_shared_operations.hpp"
#include "initialize.hpp"
#include "2d_system.hpp"

using hpx::lcos::future;
using hpx::find_here;
using hpx::lcos::wait;
using hpx::make_ready_future;
using hpx::lcos::local::dataflow;
using hpx::async;
using hpx::util::unwrap;

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;
typedef std::vector< future< shared_vec > > state_type;

typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       local_dataflow_algebra ,
                                       local_dataflow_shared_operations2d > stepper_type;

void trivial_sys( state_type &q , state_type &dpdt )
{
    for( size_t i=0 ; i<q.size() ; ++i )
    {
        dpdt[i] = dataflow( unwrap([i]( shared_vec in , shared_vec sync ) 
            { 
                hpx::cout << (boost::format("rhs %d ...\n") % i ) << hpx::flush;
                return in; 
            }) , q[i] , dpdt[i] );
    }
}

void state_swap( state_type &x1 , state_type &x2 )
{
    for( size_t i=0 ; i<x1.size() ; ++i )
    {
        future< shared_vec > tmp = x1[i];
        x1[i] = dataflow( hpx::launch::sync , 
                          unwrap([i]( shared_vec in , shared_vec sync ) { 
                                  //hpx::cout << (boost::format("swap2: %d\n") % i ) << hpx::flush;
                                  return in; }) , 
                          x2[i] , tmp );
        x2[i] = dataflow( hpx::launch::sync , 
                          unwrap([i]( shared_vec in , shared_vec sync ) { 
                                  //hpx::cout << (boost::format("swap3: %d\n") % i ) << hpx::flush;
                                  return in; }) , 
                          tmp , x1[i] );
    }
}

void state_copy( state_type &from , state_type &to )
{
    for( size_t i=0 ; i<from.size() ; ++i )
    {
        to[i] = dataflow( hpx::launch::async ,
                          unwrap( []( shared_vec from , shared_vec to ) 
            { 
                *from = *to; 
                return to;
            }),
                          from[i] , to[i] );
        from[i] = dataflow( hpx::launch::async , 
                            unwrap( []( shared_vec x , shared_vec sync ) 
            {
                return x;
            }) , 
                            from[i] , to[i] );
    }
}


int hpx_main(boost::program_options::variables_map& vm)
{

    const std::size_t N1 = vm["N1"].as<std::size_t>();
    const std::size_t N2 = vm["N2"].as<std::size_t>();
    const std::size_t G = vm["G"].as<std::size_t>();
    const bool fully_random = vm["fully_random"].as<bool>();
    const bool do_observation = vm["observe"].as<bool>();
    const std::size_t init_length = vm["init_length"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();
    const std::size_t M = N1/G;


    std::clog << "Dimension: " << N1 << "x" << N2 << ", number of rows per dataflow: " << G;
    std::clog << ", number of dataflow: " << M << ", steps: " << steps << ", dt: " << dt << std::endl;

    dvecvec p_init( N1 , dvec( N2 , 0.0 ) );

    std::uniform_real_distribution<double> distribution( -1.0 , 1.0 );
    std::mt19937 engine( 0 ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);

    if( fully_random )
    {
        for( size_t j=0 ; j<N1 ; j++ )
            std::generate( p_init[j].begin() , 
                           p_init[j].end() , 
                           std::ref(generator) );
    } else
    {
        for( size_t j=N1/2-init_length/2 ; j<N1/2+init_length/2 ; j++ )
            std::generate( p_init[j].begin()+N2/2-init_length/2 , 
                           p_init[j].begin()+N2/2+init_length/2 , 
                           std::ref(generator) );
    }

    state_type q_in( M );
    state_type p_in( M );
    state_type q_out( M );
    state_type p_out( M );

    for( size_t i=0 ; i<M ; ++i )
    {
        q_in[i] = make_ready_future( std::allocate_shared<dvecvec>( std::allocator<dvecvec>() ) );
        q_in[i] = dataflow( unwrap(initialize_zero( G , N2 )) , q_in[i] );
        q_out[i] = make_ready_future( std::allocate_shared<dvecvec>( std::allocator<dvecvec>() ) );
        q_out[i] = dataflow( unwrap(initialize_zero( G , N2 )) , q_out[i] );
        p_in[i] = make_ready_future( std::allocate_shared<dvecvec>( std::allocator<dvecvec>() ) );
        p_in[i] = dataflow( unwrap(initialize_copy( p_init , i*G , G )) , p_in[i] );
        p_out[i] = make_ready_future( std::allocate_shared<dvecvec>( std::allocator<dvecvec>() ) );
        p_out[i] = dataflow( unwrap(initialize_zero( G , N2 )) , p_out[i] );
    }

    wait( q_in );
    wait( p_in );
    std::clog.precision(10);
    std::clog << "Initialization complete, energy: " << energy( q_in , p_in ) << std::endl;

    // std::cout.precision(10);

    // //print the values
    // for( size_t m=0 ; m<M ; ++m )
    //     for( size_t g=0 ; g<G ; g++ )
    //     {
    //         for( size_t i=0 ; i<N2 ; i++ )
    //         {
    //             std::cout << (*(futures_q[m].get()))[g][i] << ',';
    //             std::cout << (*(futures_p[m].get()))[g][i] << '\t';
    //         }
    //         std::cout << std::endl;
    //     }

    hpx::util::high_resolution_timer timer;

    stepper_type stepper;

    for( size_t t=0 ; t<steps ; ++t )
    {
        auto in = std::make_pair( boost::ref(q_in) , boost::ref(p_in) );
        auto out = std::make_pair( boost::ref(q_out) , boost::ref(p_out) );
        stepper.do_step( system_2d , 
                         //trivial_sys ,
                         in ,
                         t*dt , 
                         out , 
                         dt );

        // if( do_observation && ((t%10) == 0) )
        // {
        //     obs( q_out , p_out , t*dt );
        // }

        state_swap( q_in , q_out );
        state_swap( p_in , p_out );

        // state_copy( q_out , q_in );
        // state_copy( p_out , p_in );

    }

    hpx::cout << "dataflow generation ready\n" << hpx::flush;

    wait( q_in );
    wait( p_in );

    hpx::cout << (boost::format("runtime: %fs\n") %timer.elapsed()) << hpx::flush;

    std::clog << "Integration complete, energy: " << energy( q_in , p_in ) << std::endl;

    std::cout.precision(10);

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N1",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "Dimension 1 (1024)")
        ;
    desc_commandline.add_options()
        ( "N2",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "Dimension 2 (1024)")
        ;
    desc_commandline.add_options()
        ( "G",
          boost::program_options::value<std::size_t>()->default_value(128),
          "Block size (128)")
        ;
    desc_commandline.add_options()
        ( "fully_random",
          boost::program_options::value<bool>()->default_value(false),
          "Fully random initial condition (false)")
        ;
    desc_commandline.add_options()
        ( "init_length",
          boost::program_options::value<std::size_t>()->default_value(32),
          "Initial excitation length (32)")
        ;
    desc_commandline.add_options()
        ( "observe",
          boost::program_options::value<bool>()->default_value(false),
          "Spreading observation (false)")
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
