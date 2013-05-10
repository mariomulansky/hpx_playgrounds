// Copyright 2013 Mario Mulansky
//
// some actions required to use hpx with odeint
#ifndef HPX_ODEINT_ACTIONS_HPP
#define HPX_ODEINT_ACTIONS_HPP

#include <memory>
#include <vector>
#include <algorithm>
#include <random>

typedef std::vector< std::vector<double> > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vector;

shared_vector identity_vec( shared_vector x )
{
    return x;
}

HPX_PLAIN_ACTION(identity_vec, identity_vec_action);


shared_vector initialize_2d( shared_vector x , 
                             const int size1 , const int size2 , 
                             const double value )
{
    x->resize( size1 );
    std::fill( x->begin() , x->end() , std::vector<double>( size2 , value ) );
    return x;
}

HPX_PLAIN_ACTION( initialize_2d , initialize_2d_action );

shared_vector initialize_2d_random( shared_vector x , 
                                    const int size1 , 
                                    const int size2 ,
                                    const int seed=0 )
{
    x->resize( size1 );
    std::uniform_real_distribution<double> distribution(0.0);
    std::mt19937 engine( seed ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);
    for( int i=0; i<size1 ; ++i )
    {
        (*x)[i].resize( size2 );
        std::generate( (*x)[i].begin() , (*x)[i].end() , generator );
    }
    return x;
}

HPX_PLAIN_ACTION( initialize_2d_random , initialize_2d_random_action );

shared_vector initialize_2d_from_data( shared_vector x , 
                                       const dvecvec &data ,
                                       int start ,
                                       int size
                                       )
{
    x->resize( size );
    for( int i=0; i < size ; ++i )
    {
        (*x)[i].resize( data[i+start].size() );
        std::copy( data[i+start].begin() , data[i+start].end() , (*x)[i].begin() );
    }
    return x;
}

HPX_PLAIN_ACTION( initialize_2d_from_data , initialize_2d_from_data_action );


// shared_vector initialize_random_part( shared_vector x , const int size , const int start , const int end , const int seed=0 , const double value=0.0 )
// {
//     x->resize( size );

//     std::uniform_real_distribution<double> distribution(0.0, value);
//     std::mt19937 engine( seed ); // Mersenne twister MT19937
//     auto generator = std::bind(distribution, engine);
//     //std::generate( x->begin()+start , x->begin()+end , generator );
//     return x;
// }

// HPX_PLAIN_ACTION( initialize_random_part , initialize_random_part_action );


shared_vector hpx_resize_2d( shared_vector x , shared_vector y )
{
    x->resize( y->size() );
    for( size_t i=0 ; i<x->size() ; ++i )
        (*x)[i].resize( (*y)[i].size() );
    return x;
}

HPX_PLAIN_ACTION( hpx_resize_2d , hpx_resize_2d_action );

template< typename T1 , typename T2 >
T1 sync1( T1 x , T2 sync )
{
    return x;
}

template< typename T1 , typename T2 >
struct sync1_action
    : hpx::actions::make_action<
    T1 (*)( T1 , T2 ) , 
      &sync1<T1,T2>, 
      sync1_action<T1,T2> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename T1 , typename T2 >),
    (sync1_action< T1 , T2 >))


#endif
