// Copyright 2013 Mario Mulansky
//
#include <iostream>
#include <memory>
#include <vector>

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrap.hpp>

typedef std::vector< double > dvec;
typedef std::shared_ptr<dvec> shared_vec;
typedef hpx::lcos::future< shared_vec > future_type;

struct my_init
{
    const int N;
    const double value;

    my_init( const int N_ , const double v ) 
        : N( N_ ) , value( v )
    { }

    shared_vec operator() ( shared_vec d ) const
    {
        d->resize(N);
        std::fill( d->begin() , d->end() , value );
        return d;
    }
};

int main()
{
    future_type f = hpx::make_ready_future(std::allocate_shared<dvec>(std::allocator<dvec>()));
    f = hpx::lcos::local::dataflow( hpx::util::unwrap( my_init( 1000 , 1.0 ) ) , f );

    std::cout << (*(f.get()))[0] << std::endl;

    return 0;
}
