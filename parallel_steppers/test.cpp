// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/numeric/odeint/stepper/parallel_extrapolation_stepper.hpp>
#include <boost/foreach.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp> 

#include "dataflow_wrapper.hpp"
#include "dataflow_algebra.hpp"
#include "dataflow_operations.hpp"
#include "dataflow_copy.hpp"
#include "identity_action.hpp"
//#include "serialization.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

using boost::numeric::odeint::parallel_extrapolation_stepper;

typedef dataflow_base< double > state_type;

typedef parallel_extrapolation_stepper< 8 , 
                                        state_type , double , state_type , double , 
                                        dataflow_algebra , dataflow_operations 
                                        > stepper_type;

double rhs_func( double x )
{
    hpx::cout << "rhs func waiting...\n" << hpx::flush;
    boost::this_thread::sleep( boost::posix_time::milliseconds(1000) );
    return sin(x);
}
HPX_PLAIN_ACTION( rhs_func , rhs_func_action )

void rhs( const state_type &x , state_type &dxdt , double t )
{
    dxdt = dataflow< rhs_func_action >( find_here() , x );
}

int main()
{
    state_type x_in = dataflow< identity_action >( find_here() , 1.0 );
    state_type x_out = dataflow< identity_action >( find_here() , 1.0 );
    stepper_type stepper;
    // for some reason this is necessary
    wait( x_in.get_future() );
    wait( x_out.get_future() );
    stepper.do_step( rhs , x_in , 0.0 , x_out , 0.1 );
    hpx::cout << (boost::format("%f\n") % (x_out.get_future().get())) << hpx::flush;
    return 0;
}
