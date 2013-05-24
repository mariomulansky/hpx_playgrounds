// Copyright 2013 Mario Mulansky
// resizing functionality for odeint
#ifndef HPX_DATAFLOW_WRAPPER_HPP
#define HPX_DATAFLOW_WRAPPER_HPP

#include <iostream>
#include <memory>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>

#include <hpx/include/actions.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/include/iostreams.hpp>

#include "identity_action.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

namespace boost {
namespace numeric {
namespace odeint {

template< typename T >
struct is_resizeable< dataflow_base<T> >
{
    typedef boost::false_type type;
    const static bool value = type::value;
};

template< typename T >
struct state_wrapper< dataflow_base< std::shared_ptr<T> > >
{
    typedef state_wrapper< dataflow_base<T> > state_wrapper_type;

    state_wrapper()
    {
        typedef initialize_shared_action< T > init_action;
        m_v = dataflow< init_action >( find_here() , T() );
    }

    dataflow_base< std::shared_ptr<T> > m_v;
};

} } }

#endif
