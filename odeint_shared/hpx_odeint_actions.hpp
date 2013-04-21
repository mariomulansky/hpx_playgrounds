// Copyright 2013 Mario Mulansky
//
// some actions required to use hpx with odeint
#ifndef HPX_ODEINT_ACTIONS_HPP
#define HPX_ODEINT_ACTIONS_HPP

#include <memory>

typedef std::shared_ptr< std::vector< double > > shared_vector;

shared_vector identity_vec( shared_vector x )
{
    return x;
}

HPX_PLAIN_ACTION(identity_vec, identity_vec_action);


shared_vector initialize( shared_vector x , const int size , 
                          const double value )
{
    x->resize( size );
    BOOST_FOREACH( double &elem , *x )
    {
        elem = value;
    }
    return x;
}

HPX_PLAIN_ACTION( initialize , initialize_action );


shared_vector hpx_resize( shared_vector x , const int size )
{
    x->resize( size );
    return x;
}

HPX_PLAIN_ACTION( hpx_resize , hpx_resize_action );

#endif
