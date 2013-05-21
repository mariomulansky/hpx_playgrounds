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

struct operations
{
    template< typename Value >
    struct mul 
    {
        const Value a;

        mul( const Value alpha ) : a( alpha )
        { }

        template< typename T1 , typename T2 >
        T1 operator() ( T1 x1 , T2 x2 ) const
        {
            return x1*x2*a;
        }
    };
};

int main()
{
    future_type f = hpx::make_ready_future( 1.0 );

    future_type f2 = dataflow( operations::mul<double>( 0.5 ) , f , f );

    hpx::cout << boost::format("%d\n") % f2.get() << hpx::flush;

    return 0;
}
