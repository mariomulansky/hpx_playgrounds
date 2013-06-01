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
#include <hpx/util/unwrap.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/thread/thread.hpp>

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::local::dataflow;
using hpx::find_here;
using hpx::lcos::wait;

typedef future< double > future_type;

struct operations
{
    template< typename Value >
    struct mul 
    {
        const Value a;

        mul( const Value alpha ) : a( alpha )
        { }

        double operator() ( double x1 , double x2 ) const
        {
            return x1*x2*a;
        }
    };
};


struct func
{
    double operator()( double x ) const
    {
        boost::this_thread::sleep( boost::posix_time::milliseconds(100) );
        hpx::cout << boost::format( "%f\n" ) % x << hpx::flush;
        return x*0.5;
    }
};

int main()
{
    future_type f = hpx::make_ready_future( 1.0 );

    future_type f2 = dataflow( hpx::util::unwrap(operations::mul<double>( 0.5 )) , f , f );

    hpx::cout << boost::format("%d\n") % f2.get() << hpx::flush;

    for( int n=0 ; n<10 ; ++n )
    {
        hpx::cout << boost::format("step %d\n") % n << hpx::flush;
        f = dataflow( hpx::util::unwrap( func() ) , f );
    }

    hpx::cout << boost::format("final result: %d\n") % f.get() << hpx::flush;
    
    return 0;
}
