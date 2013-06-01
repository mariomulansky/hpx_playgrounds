// Copyright 2013 Mario Mulansky
//
#include <iostream>

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrap.hpp>

using hpx::lcos::local::dataflow;

typedef hpx::lcos::future< double > future_type;

template< typename Value >
struct mul 
{
    const Value a;

    mul( const Value alpha ) : a( alpha )
    { }

    double operator() ( double x1 , double x2 ) const // this has to be const?!
    {
        return x1*x2*a;
    }
};

int main()
{
    auto functor = hpx::util::unwrap(mul<double>( 0.5 ));
    future_type f1 = hpx::make_ready_future( 1.0 );
    // compiles fine when specifically stating namespace
    //future_type f2 = hpx::lcos::local::dataflow( functor , f1 , f1 );    
    // compile fails when using dataflow without explicit namespaces
    future_type f2 = dataflow( functor , f1 , f1 );
    future_type f3 = hpx::lcos::local::dataflow( hpx::util::unwrap(mul<double>( 2.0 )) , f1 , f1 );
    
    return 0;
}
