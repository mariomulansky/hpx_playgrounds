// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/numeric/odeint.hpp>

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef dataflow_base< shared_vec > shared_dataflow;
typedef std::vector< shared_dataflow > shared_dataflow_vec;

namespace boost { namespace serialization {

template<class Archive>
void serialize(Archive & ar, shared_vec &v , const unsigned int version)
{
    // empty for now
}

} }

shared_vec initialize( shared_vec x , const int size , 
                       const double value )
{
    std::cout << "initializing..." << std::endl;
    x->resize( size );
    std::cout << "resizing complete" << std::endl;
    BOOST_FOREACH( double &elem , *x )
    {
        elem = value;
    }
    std::cout << "initialization complete" << std::endl;
    return x;
}

HPX_PLAIN_ACTION( initialize , initialize_action );

int main()
{
    shared_dataflow x = dataflow<initialize_action>( find_here() , 
                                                     std::allocate_shared<dvec>( std::allocator<dvec>() ) , 
                                                     10 , 0.0 );
    x.get_future().wait();

    const int N = 10;

    shared_dataflow_vec data( N );
    std::vector< future<shared_vec> > futures( N );
    for( size_t i=0 ; i<N ; ++i )
    {
        data[i] = dataflow<initialize_action>( find_here() , 
                                                     std::allocate_shared<dvec>( std::allocator<dvec>() ) , 
                                                     10 , 0.0 );
        futures[i] = data[i].get_future();
    }

    wait( futures );

    return 0;
}
