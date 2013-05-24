#ifndef DATAFLOW_COPY
#define DATAFLOW_COPY

#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/include/iostreams.hpp>

#include "identity_action.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

namespace boost { namespace numeric { namespace odeint {

template< size_t N>
struct copy_impl< dataflow_base< std::shared_ptr< std::array<double,N> > >, 
                  dataflow_base< std::shared_ptr< std::array<double,N> > > >
{
    typedef dataflow_base< std::shared_ptr< std::array<double,N> > > state_type;

    static void copy( const state_type &from , state_type &to )
    {
        typedef copy_shared_container_action< std::array<double,N> > copy_action;
        to = dataflow< copy_action >( find_here() , from , to );
    }
};

} } }

#endif
