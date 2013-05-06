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

#include "dataflow_shared_resize.hpp"
#include "dataflow_shared_algebra.hpp"
#include "dataflow_shared_operations.hpp"
#include "serialization.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

using hpx::cout;
using hpx::flush;

using boost::numeric::odeint::euler;
using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;

typedef dataflow_base< double > df_base;
typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef std::vector< dataflow_base< shared_vec > > state_type;
typedef std::vector< dvec > vec_dvec;

const double KAPPA=3.5;
const double LAMBDA=4.5;

//typedef euler< state_type , double , state_type , double , dataflow_shared_algebra , dataflow_operations > stepper_type;
//typedef euler< state_type > stepper_type;
typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       dataflow_shared_algebra ,
                                       dataflow_operations > stepper_type;

shared_vec osc_operation( shared_vec q , shared_vec q_l , shared_vec q_r , shared_vec dpdt )
{
    using std::pow;
    using std::abs;
    using boost::math::sign;

    const size_t M = q->size();

    double coupling = pow( abs((*q_l)[q_l->size()-1]-(*q)[0]) , LAMBDA-1 ) * sign((*q_l)[q_l->size()-1]-(*q)[0]);

    // (*dpdt)[0] = -pow( abs((*q)[0]) , KAPPA-1 ) * sign( (*q)[0] )
    //     - pow( abs((*q)[0]-(*q)[1]) , LAMBDA-1 ) * sign( (*q)[0]-(*q)[1] )
    //     - pow( abs((*q)[0]-(*q_l)[q_l->size()-1]) , LAMBDA-1 ) * sign((*q)[0]-(*q_l)[q_l->size()-1]);

    for( size_t i=0 ; i<M-1 ; ++i )
    {
        (*dpdt)[i] = -pow( abs((*q)[i]) , KAPPA-1 ) * sign((*q)[i])
            + coupling;
        coupling = pow( abs((*q)[i]-(*q)[i+1]) , LAMBDA-1 ) * sign((*q)[i]-(*q)[i+1]);
        (*dpdt)[i] -= coupling;
    }

    (*dpdt)[M-1] = -pow( abs((*q)[M-1]) , KAPPA-1 ) * sign((*q)[M-1]) 
        - pow( abs((*q)[M-1] - (*q_r)[0]) , LAMBDA-1 ) * sign((*q)[M-1]-(*q_r)[0])
        + coupling;

    return dpdt;
}

HPX_PLAIN_ACTION(osc_operation, osc_operation_action);

void osc_sys( const state_type &q , state_type &dpdt )
{
    const size_t N = q.size();
    //std::cout << "system call with N=" << N << std::endl;
    dpdt[0] = dataflow< osc_operation_action >( find_here() ,
                                                q[0] , q[N-1] , q[1] ,
                                                dpdt[0] );

    for( size_t i=1 ; i<N-1 ; ++i )
        dpdt[i] = dataflow< osc_operation_action >( find_here() , 
                                                    q[i] , q[i-1] , q[i+1] ,
                                                    dpdt[i] );

    dpdt[N-1] = dataflow< osc_operation_action >( find_here() ,
                                                  q[N-1] , q[N-2] , q[0] ,
                                                  dpdt[N-1] );
    //std::cout << "-------" << std::endl;
}


shared_vec simple_osc_operation( shared_vec q , shared_vec dpdt )
{
    using boost::math::pow;

    const size_t M = q->size();

    //std::cout << M << std::endl;

    (*dpdt)[0] = -pow<3>((*q)[0]);

    for( size_t i=1 ; i<M-1 ; ++i )
    {
        (*dpdt)[i] = -pow<3>((*q)[i]);
    }

    (*dpdt)[M-1] = -pow<3>((*q)[M-1]);
    
    return dpdt;
}

HPX_PLAIN_ACTION(simple_osc_operation, simple_osc_operation_action);

void simple_osc_sys( const state_type &q , state_type &dpdt )
{
    const size_t N = q.size();
    //std::cout << "system call with N=" << N << std::endl;
    dpdt[0] = dataflow< simple_osc_operation_action >( find_here() ,
                                                q[0] ,
                                                dpdt[0] );

    for( size_t i=1 ; i<N-1 ; ++i )
        dpdt[i] = dataflow< simple_osc_operation_action >( find_here() , 
                                                    q[i] ,
                                                    dpdt[i] );

    dpdt[N-1] = dataflow< simple_osc_operation_action >( find_here() ,
                                                  q[N-1] ,
                                                  dpdt[N-1] );
    //std::cout << "-------" << std::endl;
}


shared_vec sync_identity( shared_vec x , shared_vec sync1 , shared_vec sync2 , shared_vec sync3 )
{
    return x;
}

HPX_PLAIN_ACTION( sync_identity , sync_identity_action )

shared_vec sync_identity2( shared_vec x , shared_vec sync1 )
{
    return x;
}

HPX_PLAIN_ACTION( sync_identity2 , sync_identity2_action )

// syncronized swap emulating nearest neighbor coupling
void synchronized_swap( state_type &x_in , state_type &x_out )
{
    const size_t N = x_in.size();
    dataflow_base< shared_vec > x_0 = x_out[0];
    dataflow_base< shared_vec > x_left = x_out[0];
    dataflow_base< shared_vec > x_tmp = dataflow< sync_identity_action >( find_here() , x_in[0] , 
                                                                          x_out[0] , x_out[N-1] , x_out[1] );
    x_in[0] = dataflow< sync_identity2_action >( find_here() , x_out[0] , x_tmp );
    x_out[0] = x_tmp;
    for( size_t n=1 ; n<N-1 ; ++n )
    {
        x_tmp = dataflow< sync_identity_action >( find_here() , x_in[n] , 
                                              x_out[n] , x_left , x_out[n+1] );
        x_left = x_out[n];
        x_in[n] = dataflow< sync_identity2_action >( find_here() , x_out[n] , x_tmp );
        x_out[n] = x_tmp;
    }
    x_tmp = dataflow< sync_identity_action >( find_here() , x_in[N-1] , 
                                          x_out[N-1] , x_left , x_0 );
    x_in[N-1] = dataflow< sync_identity2_action >( find_here() , x_out[N-1] , x_tmp );
    x_out[N-1] = x_tmp;
}

double energy( const dvec &q , const dvec &p )
{
    using std::pow;
    using std::abs;
    using boost::math::pow;
    const size_t N=q.size();
    double e(0.0);
    for( size_t n=0 ; n<N-1 ; ++n )
    {
        e += 0.5*pow<2>( p[n] ) 
            + pow( abs(q[n]) , KAPPA )/KAPPA 
            + pow( abs(q[n]-q[n+1]) , LAMBDA )/LAMBDA;
    }
    e += 0.5*pow<2>( p[N-1] ) 
        + pow( abs(q[N-1]) , KAPPA )/KAPPA 
        + pow( abs(q[N-1]-q[0]) , LAMBDA )/LAMBDA;
    return e;
}

double simple_energy( const dvec &q , const dvec &p )
{
    using boost::math::pow;
    const size_t N=q.size();
    double e(0.0);
    for( size_t n=0 ; n<N ; ++n )
    {
        e += 0.5*pow<2>( p[0] ) 
            + pow<4>( q[0] )/4.0;
    }
    return e;
}

double energy( const std::vector< future<shared_vec> > &q_futures ,
               const std::vector< future<shared_vec> > &p_futures )
{
    const size_t N=q_futures.size();

    dvec q , p;
    
    for( size_t i=0 ; i<N ; ++i )
    {
        shared_vec data = q_futures[i].get();
        q.insert( q.end() , data->begin() , data->end() );
        data = p_futures[i].get();
        p.insert( p.end() , data->begin() , data->end() );
    }
    return energy( q , p );
}

int hpx_main(boost::program_options::variables_map& vm)
{

    const std::size_t N_ = vm["N"].as<std::size_t>();
    const std::size_t M = vm["granularity"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();

    const std::size_t N = N_/M;

    // std::clog << "number of dataflows: " << N << ", number of elements per dataflow: " << M;
    // std::clog <<  ", steps: " << steps << ", dt: " << dt << std::endl;

    double mean_t( 0.0 );

    for( size_t n=0 ; n<10 ; ++n )
    {

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
        // std::clog << "Initialization complete, energy: " << energy( futures_q , futures_p ) << std::endl;

        hpx::util::high_resolution_timer timer;
    
        // do the numerical integration using dataflow objects
        stepper_type stepper;
        for( size_t t=0 ; t<steps ; ++t )
        {
            auto in = std::make_pair( boost::ref(q_in) , boost::ref(p_in) );
            auto out = std::make_pair( boost::ref(q_out) , boost::ref(p_out) );
            stepper.do_step( osc_sys , 
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

        mean_t += timer.elapsed();
    }

    cout << (boost::format("%d") % (N_/M) ) << '\t' << (boost::format("%d") % M);
    cout << '\t' << (boost::format("%f") % (mean_t/10)) << "\n" << flush;
    // std::clog << "Final energy: " << energy( futures_q , futures_p ) << std::endl;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "Number of elements (1024)")
        ;
    desc_commandline.add_options()
        ( "granularity",
          boost::program_options::value<std::size_t>()->default_value(128),
          "Granularity - Number of elements in each dataflow (100)")
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
