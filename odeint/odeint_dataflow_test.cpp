// Copyright 2013 Mario Mulansky

#include <iostream>
#include <vector>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/numeric/odeint.hpp>

#include "dataflow_algebra.hpp"
#include "dataflow_operations.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

using boost::numeric::odeint::euler;

typedef dataflow_base< double > df_base;
typedef std::vector< df_base > state_type;

typedef euler< state_type , double , state_type , double , dataflow_algebra , dataflow_operations > stepper_type;

const int N = 1000;
const size_t steps = 100;
const double dt = 0.01;
const double lmbd = 0.01;

double rhs_operation( const double x )
{
    return -lmbd*x;
}

HPX_PLAIN_ACTION(rhs_operation, rhs_operation_action);

double identity(double initial_value)
{
    return initial_value;
}

HPX_PLAIN_ACTION(identity, identity_action);

void sys( const state_type &x , state_type &dxdt , const double dt )
{
    const size_t N = x.size();
    for( size_t i=0 ; i<N ; ++i )
        dxdt[i] = dataflow< rhs_operation_action >( find_here() , x[i] );
}

int main()
{
    // boost::program_options::options_description
    //    desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    // desc_commandline.add_options()
    //     ( "N",
    //       boost::program_options::value<std::size_t>()->default_value(100),
    //       "System size")
    //     ;
    // desc_commandline.add_options()
    //     ( "steps",
    //       boost::program_options::value<std::size_t>()->default_value(100),
    //       "steps")
    //     ;

    state_type x( N , dataflow< identity_action >( find_here() , 100.0 ) );
    //df_base t = dataflow< identity_action >( find_here() , 0.0 );
    //df_base dt = dataflow< identity_action >( find_here() , 0.1 );
    double t = 0.0;

    stepper_type stepper;

    hpx::util::high_resolution_timer timer;

    boost::numeric::odeint::integrate_n_steps( stepper , sys , x , t , dt , steps );

    std::clog << "Dependency tree built" << std::endl;

    std::vector< future<double> > futures( N );
    for( int i=0 ; i<N ; ++i )
        futures[i] = x[i].get_future();
    wait( futures );

    std::cout << timer.elapsed() << std::endl;
    return 0;
}
