// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/numeric/odeint.hpp>

#include "hpx_dataflow_resize.hpp"
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
typedef std::vector< double > dvec;
typedef std::vector< dataflow_base< dvec > > state_type;

typedef euler< state_type , double , state_type , double , dataflow_sub_algebra , dataflow_operations > stepper_type;

const double lmbd = 0.01;

dvec rhs_operation( dvec x , dvec dxdt )
{
    for( size_t i=0 ; i<x.size() ; ++i )
        dxdt[i] = -x[i]; //pow(abs(x[i]),0.5) * sin(x[i]) + cos(x[i]);
    return dxdt;
}

HPX_PLAIN_ACTION(rhs_operation, rhs_operation_action);

void sys( const state_type &x , state_type &dxdt , const double dt )
{
    const size_t N = x.size();
    for( size_t i=0 ; i<N ; ++i )
        dxdt[i] = dataflow< rhs_operation_action >( find_here() , x[i] , dxdt[i] );
}

int hpx_main(boost::program_options::variables_map& vm)
{

    const std::size_t N = vm["N"].as<std::size_t>();
    const std::size_t M = vm["M"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();

    std::clog << "number of dataflows: " << N << ", number of elements per dataflow: " << M;
    std::clog <<  ", steps: " << steps << ", dt: " << dt << std::endl;

    state_type x( N );
    for( size_t i=0 ; i<N ; ++i )
        x[i] = dataflow< identity_vec_action >( find_here() , dvec( M , 100.0 ) );
    double t = 0.0;

    hpx::util::high_resolution_timer timer;
    
    // do the numerical integration using dataflow objects
    boost::numeric::odeint::integrate_n_steps( stepper_type() , sys , x , t , dt , steps );

    std::clog << "Dependency tree built" << std::endl;

    std::vector< future<dvec> > futures( N );
    for( size_t i=0 ; i<N ; ++i )
    {
        futures[i] = x[i].get_future();
    }
    
    // here we wait for the results
    wait( futures );

    std::clog << "Calculation finished in " << timer.elapsed() << "s" << std::endl;

    // print the values
    // for( size_t i=0 ; i<N ; ++i )
    //     for( size_t m=0 ; m<M ; m++ )
    //         std::cout << futures[i].get()[m] << '\t';
    // std::cout << std::endl;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N",
          boost::program_options::value<std::size_t>()->default_value(4),
          "Number of dataflows (4)")
        ;
    desc_commandline.add_options()
        ( "M",
          boost::program_options::value<std::size_t>()->default_value(10000),
          "Number of elements in each dataflow (100)")
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
