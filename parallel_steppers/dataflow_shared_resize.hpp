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

namespace boost {
namespace numeric {
namespace odeint {

template< typename T , size_t N >
struct is_resizeable< boost::array<T,N> >
{
    typedef boost::false_type type;
    const static bool value = type::value;
};

} } }


#endif
