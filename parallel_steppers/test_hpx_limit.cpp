// Copyright Mario Mulansky 2013

#include <iostream>

// anything > 10 here gives boost MPL errors
#define HPX_LIMIT 10

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

//define large action
double func( double x1 , double x2 , double x3 , double x4 , double x5 , double x6 , double x7 )
{
    return x1;
}

HPX_PLAIN_ACTION( func , large_action );

int main()
{
    // action too big to compile...
    dataflow_base< double > df = dataflow< large_action >( find_here() , 
                                                           1 , 1 , 1 , 1 , 1 , 1 , 1 );
    std::cout << df.get_future().get() << std::endl;
    return 0;
};
