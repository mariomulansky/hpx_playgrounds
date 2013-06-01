#ifndef IDENTITY_ACTION_HPP
#define IDENTITY_ACTION_HPP

#include <memory>

#include <hpx/hpx.hpp>
#include <hpx/include/actions.hpp>

double identity( const double x )
{
    //hpx::cout << "identity action\n" << hpx::flush;
    return x;
}
HPX_PLAIN_DIRECT_ACTION( identity , identity_action );

template< typename T>
T identity_tmpl( const T &x )
{
    hpx::cout << "identity action\n" << hpx::flush;
    return x;
}
template< typename T >
struct identity_tmpl_action
    : hpx::actions::make_direct_action<
    T (*)( const T& ) , 
    &identity_tmpl<T>, 
    identity_tmpl_action<T> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename T >),
    (identity_tmpl_action< T >))


template< typename T>
std::shared_ptr<T> initialize_shared( const T &x )
{
    //hpx::cout << "identity action\n" << hpx::flush;
    std::shared_ptr<T> p = std::allocate_shared< T >( std::allocator<T>() );
    std::copy( x.begin() , x.end() , p->begin() );
    return p;
}
template< typename T >
struct initialize_shared_action
    : hpx::actions::make_direct_action<
    std::shared_ptr<T> (*)( const T& ) , 
    &initialize_shared<T>, 
    initialize_shared_action<T> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename T >),
    (initialize_shared_action< T >))


template< typename T1 , typename T2 >
T1 identity_sync1( const T1 &x , const T2 &sync )
{
    return x;
}

template< typename T1 , typename T2 >
struct identity_sync1_action
    : hpx::actions::make_direct_action<
    T1 (*)( const T1& , const T2 & ) , 
        &identity_sync1<T1,T2>, 
        identity_sync1_action<T1,T2> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
(template< typename T1 , typename T2 >),
(identity_sync1_action< T1 , T2 >))


template< typename T1 , typename T2 >
T1 identity_sync2( const T1 &x , const T2 &sync1 , const T2 &sync2 )
{
    //hpx::cout << (boost::format( "%d\t%d\t%d\n" ) % (sync1->size()) % (sync2->size()) % (sync3->size()) ) << hpx::flush;
    return x;
}

template< typename T1 , typename T2 >
struct identity_sync2_action
    : hpx::actions::make_direct_action<
    T1 (*)( const T1& , const T2& , const T2& ) ,
        &identity_sync2<T1,T2>, 
        identity_sync2_action<T1,T2> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
(template< typename T1 , typename T2 >),
(identity_sync2_action< T1 , T2 >))


template< typename T1 , typename T2 >
T1 identity_sync3( const T1 &x , const T2 &sync1 , const T2 &sync2 , const T2 &sync3 )
{
    //hpx::cout << (boost::format( "%d\t%d\t%d\n" ) % (sync1->size()) % (sync2->size()) % (sync3->size()) ) << hpx::flush;
    return x;
}

template< typename T1 , typename T2 >
struct identity_sync3_action
    : hpx::actions::make_direct_action<
    T1 (*)( const T1& , const T2& , const T2& , const T2& ) , 
        &identity_sync3<T1,T2>, 
        identity_sync3_action<T1,T2> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
(template< typename T1 , typename T2 >),
(identity_sync3_action< T1 , T2 >))


template< typename Container >
std::shared_ptr< Container > 
copy_shared_container( const std::shared_ptr< Container > from , 
                       std::shared_ptr< Container > to )
{
    std::copy( from->begin() , from->end() , to->begin() );
    return to;
}

template< typename Container >
struct copy_shared_container_action
    : hpx::actions::make_direct_action<
    std::shared_ptr<Container> (*)( const std::shared_ptr< Container > , 
                                    std::shared_ptr< Container > ) , 
        &copy_shared_container<Container>, 
        copy_shared_container_action<Container> >
{};

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
(template< typename Container >),
(copy_shared_container_action< Container >))

#endif
