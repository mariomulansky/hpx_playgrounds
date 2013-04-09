// Copyright 2013 Mario Mulansky

#include <iostream>
#include <vector>

// When using the dataflow component we have to define the following constant
// as this component uses up to 6 arguments for one of its components.
#define HPX_LIMIT 6

#include <hpx/hpx_main.hpp>
#include <hpx/components/dataflow/dataflow.hpp>

#include "euler.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;

typedef dataflow_base<double> df_base;
typedef std::vector< df_base > state_type;

double osc_op1( double x )
{
    return x;
}

double osc_op2( double x )
{
    return -sin(x);
}

HPX_PLAIN_ACTION(osc_op1, osc_op1_action);
HPX_PLAIN_ACTION(osc_op2, osc_op2_action);

double identity(double initial_value)
{
    return initial_value;
}

HPX_PLAIN_ACTION(identity, identity_action);


void osc( state_type &x , state_type &dxdt , hpx::naming::id_type loc )
{
    dxdt[0] = dataflow< osc_op1_action >( loc , x[1] );
    dxdt[1] = dataflow< osc_op2_action >( loc , x[0] );
}

const double dt = 0.01;
const size_t steps = 1000;

int main()
{
    hpx::naming::id_type here = hpx::find_here();
    
    state_type x( 2 );
    x[0] = dataflow< identity_action >( here , 1.0 );
    x[1] = dataflow< identity_action >( here , 0.0 );
    state_type dxdt( 2 );
  
    for( size_t t=0 ; t<steps ; ++t )
    {
        osc( x , dxdt , here );
        hpx_step::euler_step( x , dxdt , here );
        if( steps % 10 == 0 )
            std::cout << x[0].get_future().get() << '\t' << x[1].get_future().get() << std::endl;
    }
    
    return 0;
};
