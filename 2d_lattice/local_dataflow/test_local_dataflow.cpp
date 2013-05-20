// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::local::dataflow;
using hpx::find_here;
using hpx::lcos::wait;

typedef future< double > future_type;

struct mul 
{
    const double a;

    mul( const double alpha ) : a( alpha )
    { }

    double operator() ( const double x )
    {
        return x*a;
    }
};

int main()
{
    
    future_type f = hpx::make_ready_future( 1.0 );

    future_type f2 = dataflow( mul( 0.5 ) , f );

    hpx::cout << boost::format("%d\n") % f2.get() << hpx::flush;

    return 0;
}
