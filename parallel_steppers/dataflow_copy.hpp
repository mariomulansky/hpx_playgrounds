#ifndef DATAFLOW_COPY
#define DATAFLOW_COPY

#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/include/iostreams.hpp>

#include "identity_action.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

namespace boost { namespace numeric { namespace odeint {

template< typename T >
struct copy_impl< dataflow_base<T> , dataflow_base<T> >
{
    static void copy( const dataflow_base<T> &from , dataflow_base<T> &to )
    {
        hpx::cout << "copy\n" << hpx::flush;
        typedef identity_tmpl_action< T > ident;
        to = dataflow< ident >( find_here() , from );
    }
};

} } }

#endif
