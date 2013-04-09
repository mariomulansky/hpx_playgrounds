#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/foreach.hpp>

#include <hpx/hpx_main.hpp>
#include <hpx/include/lcos.hpp>

#include "transform.hpp"

typedef std::vector<double> dvec;

const int N = 10000000;

struct add {

    template <typename>
    struct result;

    template <typename F, typename T1 , typename T2>
    struct result<F(T1,T2)>{
        typedef void type;
    };
    
    void operator()( double &x1 , const double x2 ) const
    {
        x1 += x2;
    }

};


struct complex_op {

    template <typename>
    struct result;

    template <typename F, typename T1 , typename T2>
    struct result<F(T1,T2)>{
        typedef void type;
    };
    
    void operator()( double &x1 , const double x2 ) const
    {
        x1 += 1.0/x2 * sin( x1 ) - pow( x2 , 3.5 ) + sqrt( x1 + 0.2 ) + 1.0/pow( x1 , 5.3 ) - cos( x2 );
    }

};

namespace sequential {
    // non-parallel version of the transform_inplace
    template< typename F >
    void transform_inplace( dvec &v1 , dvec const &v2 , F f )
    {
        for( std::size_t i=0 ; i<v1.size() ; ++i )

            {
                f( v1[i] , v2[i] );
            }
    }
}

int main()
{
    dvec v1( N ) , v2( N );
    
    for( int i=0 ; i<N ; ++i )
    {
        v1[i] = i;
        v2[i] = 1.0/(i+1.0);
        //std::cout << v1[i] << " + " << v2[i] << std::endl;
    }

    // schedules the transformation
    std::vector< lcos::future<void> > futures = hpx::transform_inplace_flat( v1 , v2 , complex_op() );

    std::clog << "Start calculation with " << futures.size() << " chunks" << std::endl;

    hpx::util::high_resolution_timer t;

    // run the transformation
    hpx::lcos::wait( futures );

    //sequential::transform( v1 , v2 , complex_op() );

    std::clog << "Result after: " << t.elapsed() << "s" << std::endl;


    std::cout << v1[N/2] << std::endl;
    // BOOST_FOREACH( double &x , v1 )
    // {
    //    std::cout << x << std::endl;
    // }

    return 0;
}
