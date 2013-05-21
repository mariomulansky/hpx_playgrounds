#include <iostream>

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/lcos/future.hpp>
//#include <hpx/lcos/local/packaged_task.hpp>
#include <hpx/lcos/local/dataflow.hpp>

using hpx::make_ready_future;
using hpx::lcos::future;
using hpx::lcos::local::dataflow;

typedef future<double> fut_doub;

void test( fut_doub &res , /*const*/ fut_doub in ) // this shall not be const, why?
{
    res = dataflow( []( const double d ) { return 2*d; } , in );
}


int main()
{
    fut_doub d = hpx::make_ready_future( 1.0 );
    fut_doub res;
    test( res , d );
    return 0;
}
