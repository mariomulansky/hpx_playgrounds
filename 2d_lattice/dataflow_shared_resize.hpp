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

#include <hpx/include/actions.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/include/util.hpp>
#include <hpx/lcos/async.hpp>

#include "hpx_odeint_actions.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

typedef dataflow_base< double > df_base;
typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;
typedef std::vector< dataflow_base< shared_vec > > state_type;

namespace boost {
namespace numeric {
namespace odeint {

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
        return ( ( x1.size() == x2.size() ) );
    }
};

template<>
struct resize_impl< state_type , state_type >
{
    static void resize( state_type &x1 ,
                        const state_type &x2 )
    {
        // allocate required memory
        x1.resize( x2.size() );
        for( size_t i=0 ; i < x2.size() ; ++i )
        {
            x1[i] = dataflow< hpx_resize_2d_action >( find_here() , 
                                                      std::allocate_shared<dvecvec>( std::allocator<dvecvec>() ) , 
                                                      x2[i].get_future().get() );
        }
    }
};

} } }


#endif
