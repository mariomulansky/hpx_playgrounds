// Copyright 2013 Mario Mulansky
// resizing functionality for odeint
#ifndef HPX_DATAFLOW_SHARED_RESIZE_HPP
#define HPX_DATAFLOW_SHARED_RESIZE_HPP

#include <iostream>
#include <memory>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resize.hpp>
#include <boost/numeric/odeint/util/same_size.hpp>

#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrap.hpp>

using hpx::lcos::future;
using hpx::make_ready_future;
using hpx::lcos::local::dataflow;
using hpx::util::unwrap;

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;
typedef std::vector< future< shared_vec > > state_type;

namespace boost {
namespace numeric {
namespace odeint {

    /*
template<>
struct state_wrapper< state_type >
{
    state_wrapper()
    {
        m_v = make_ready_future( state_type() );
    }

    state_type m_v;
};
    */

template<>
struct is_resizeable< state_type >
{
    typedef boost::true_type type;
    const static bool value = type::value;
};

template<>
struct same_size_impl< state_type , state_type >
{
    static bool same_size( const state_type &x1 ,
                           const state_type &x2 )
    {
        // not quite complete...
        return ( ( x1.size() == x2.size() ) );
    }
};

template<>
struct resize_impl< state_type , state_type >
{
    static void resize( state_type &x1 ,
                        state_type x2 )
    {
        //std::cout << "resizing..." << std::endl;
        // allocate required memory
        x1.resize( x2.size() );
        for( size_t i=0 ; i < x2.size() ; ++i )
        {
            x1[i] = dataflow( unwrap([]( shared_vec v2 )
                              {
                                  shared_vec tmp = std::allocate_shared<dvecvec>( std::allocator<dvecvec>() );
                                  tmp->resize( v2->size() );
                                  for( size_t n=0 ; n<v2->size() ; ++n )
                                      (*tmp)[n].resize( (*v2)[n].size() );
                                  return tmp;
                              }) ,
                              x2[i] );
        }
        //std::cout << "resizing complete" << std::endl;
    }
};

} } }


#endif
