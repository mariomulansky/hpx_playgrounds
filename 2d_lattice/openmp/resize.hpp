// Copyright Mario Mulansky 2013

#include <vector>
#include <iostream>

#include <boost/numeric/odeint/util/resize.hpp>

namespace boost { namespace numeric { namespace odeint {

typedef std::vector< std::vector< double > > state_type;

template<>
struct resize_impl< state_type , state_type >
{
    static void resize( state_type &out , const state_type &in )
    {
        out.resize( boost::size( in ) );
        typename state_type::iterator begin1 = boost::begin( out );
        typename state_type::const_iterator begin2 = boost::begin( in );
        while( begin1 != boost::end( out ) )
            (*begin1++).resize( boost::size( *begin2++ ) );
    }
};

} } } 
