// Copyright 2013 Mario Mulansky
//
// some actions required to use hpx with odeint
#ifndef HPX_ODEINT_ACTIONS_HPP
#define HPX_ODEINT_ACTIONS_HPP

#include <memory>
#include <vector>
#include <algorithm>
#include <random>

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
    std::fill( x->begin() , x->end() , value );
    return x;
}

HPX_PLAIN_ACTION( initialize , initialize_action );

shared_vector initialize_random( shared_vector x , const int size , const int seed=0 )
{
    x->resize( size );

    std::uniform_real_distribution<double> distribution(0.0);
    std::mt19937 engine( seed ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);
    std::generate( x->begin() , x->begin() , generator );
    return x;
}

HPX_PLAIN_ACTION( initialize_random , initialize_random_action );

shared_vector initialize_random_part( shared_vector x , const int size , const int start , const int end , const int seed=0 , const double value=0.0 )
{
    x->resize( size );

    std::uniform_real_distribution<double> distribution(0.0, value);
    std::mt19937 engine( seed ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);
    std::generate( x->begin()+start , x->begin()+end , generator );
    return x;
}

HPX_PLAIN_ACTION( initialize_random_part , initialize_random_part_action );


shared_vector hpx_resize( shared_vector x , const int size )
{
    x->resize( size );
    return x;
}

HPX_PLAIN_ACTION( hpx_resize , hpx_resize_action );

#endif
