//  Copyright (c) 2013 Mario Mulansky
// 
//  Distributed under the Boost Software License, Version 1.0. (See accompanying 
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// This test case demonstrates the issue described in #775: runtime error with 
// local dataflow (copying futures?).

#include <iostream>

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrap.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/thread/thread.hpp>

using hpx::lcos::local::dataflow;
using hpx::util::unwrap;

typedef hpx::lcos::future< double > future_type;

struct mul
{
    double operator()( double x1 , double x2 ) const
    {
        //boost::this_thread::sleep( boost::posix_time::milliseconds(1000) );
        //hpx::cout << boost::format( "mul: %f , %f\n" ) % x1 %x2 << hpx::flush;
        return x1*x2;
    }
};

struct divide
{
    double operator()( double x1 , double x2 ) const
    {
        //boost::this_thread::sleep( boost::posix_time::milliseconds(1000) );
        //hpx::cout << boost::format( "div: %f , %f\n" ) % x1 %x2 << hpx::flush;
        return x1/x2;
    }
};

void future_swap( future_type &f1 , future_type &f2 )
{
    //future_type tmp = hpx::lcos::local::dataflow( unwrap( []( double x ){ return x; } ) , f1 );
    future_type tmp = f1;
    f1 = dataflow( hpx::launch::sync ,
                   unwrap( []( double x , double sync ){ return x; } ) ,
                   f2 , f1 );
    f2 = dataflow( hpx::launch::sync ,
                   unwrap( []( double x , double sync ){ return x; } ) ,
                   tmp , f1 );
}

void future_copy( future_type &f1 , future_type &f2 , int i)
{
    f1 = dataflow( hpx::launch::async , unwrap( [i]( double x , double sync )
        { 
            hpx::cout << (boost::format("copy %d\n") % i) << hpx::flush;
            return x; 
        } ) , f2 , f1 );
}

int main()
{
    future_type f1 = hpx::make_ready_future( 2.0 );
    future_type f2 = hpx::make_ready_future( 3.0 );

    for( int n=0 ; n<20 ; ++n )
    {
        f1 = dataflow( hpx::launch::async , unwrap(mul()) , f1 , f2 );
        f2 = dataflow( hpx::launch::async , unwrap(divide()) , f1 , f2 );
        future_swap( f1 , f2 );
        //wait( f1 );
    }

    hpx::cout << "futures ready\n" << hpx::flush;

    hpx::cout << boost::format("f1: %d\n") % f1.get() << hpx::flush;
    hpx::cout << boost::format("f2: %d\n") % f2.get() << hpx::flush;
   
    return 0;
}
