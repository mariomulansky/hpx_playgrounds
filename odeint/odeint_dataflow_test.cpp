// Copyright 2013 Mario Mulansky
//
// trivial example on how to use Boost.odeint with hpx dataflows

#include <iostream>
#include <vector>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
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

int hpx_main(boost::program_options::variables_map& vm)
{

    const std::size_t N = vm["N"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();

    std::clog << "system size: " << N << ", steps: " << steps << ", dt: " << dt << std::endl;

    state_type x( N , dataflow< identity_action >( find_here() , 100.0 ) );
    double t = 0.0;

    hpx::util::high_resolution_timer timer;
    
    // do the numerical integration using dataflow objects
    boost::numeric::odeint::integrate_n_steps( stepper_type() , sys , x , t , dt , steps );

    std::clog << "Dependency tree built" << std::endl;

    std::vector< future<double> > futures( N );
    for( size_t i=0 ; i<N ; ++i )
    {
        futures[i] = x[i].get_future();
    }
    
    // here we wait for the results
    wait( futures );

    std::clog << "Calculation finished in " << timer.elapsed() << "s" << std::endl;

    // print the values
    // for( size_t i=0 ; i<N ; ++i )
    //     std::cout << futures[i].get() << '\t';
    // std::cout << std::endl;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N",
          boost::program_options::value<std::size_t>()->default_value(100),
          "System size")
        ;
    desc_commandline.add_options()
        ( "steps",
          boost::program_options::value<std::size_t>()->default_value(100),
          "time steps")
        ;
    desc_commandline.add_options()
        ( "dt",
          boost::program_options::value<double>()->default_value(0.01),
          "step size")
        ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
