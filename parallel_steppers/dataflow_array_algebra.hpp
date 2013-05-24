// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_ARRAY_ALGEBRA_HPP
#define DATAFLOW_ARRAY_ALGEBRA_HPP

#include <hpx/include/iostreams.hpp>

#include "dataflow_operations.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::lcos::future;
using hpx::find_here;

namespace hpx_odeint_actions {

    template< typename S , typename Operation >
    S operation_3( S x1 , const S &x2 , const S &x3 , 
                   Operation op )
    {
        hpx::cout << (boost::format("size: %d , %d , %d \n") % (x1.size()) % (x2.size()) % (x3.size()) ) << hpx::flush;
        for( size_t i=0 ; i<x1.size() ; ++i )
            op( x1[i] , x2[i] , x3[i] );
        return x1;
    }


    template< typename S , typename Operation >
    struct operation_3_action
        : hpx::actions::make_action<
        S (*)( S , 
               const S& , 
               const S& , 
               Operation op ) , 
        &operation_3<S,Operation>, 
        operation_3_action<S,Operation> >
    //        boost::mpl::true_ >
    {};


    template< typename S , typename Operation >
    S operation_4( S x1 , const S &x2 , const S &x3 , const S &x4 ,
                   Operation op )
    {
        //hpx::cout << (boost::format("vals: %d , %d , %d , %d \n") % (x1) % (x2) % (x3) , (x4) ) << hpx::flush;
        for( size_t i=0 ; i<x1.size() ; ++i )
            op( x1[i] , x2[i] , x3[i] , x4[i] );
        return x1;
    }


    template< typename S , typename Operation >
    struct operation_4_action
        : hpx::actions::make_action<
        S (*)( S , 
               const S& , 
               const S& ,
               const S& ,
               Operation op ) , 
        &operation_4<S,Operation>, 
        operation_4_action<S,Operation> >
    //        boost::mpl::true_ >
    {};


    template< typename S >
    S identity_sync( S x , const S &sync )
    {
        //hpx::cout << "sync\n" << hpx::flush;
        return x;
    }

    template< typename S >
    struct identity_sync_action
        : hpx::actions::make_direct_action<
        S (*)( S , 
               const S& ) , 
        &identity_sync<S>, 
        identity_sync_action<S> >
    {};

}

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename S , typename Operation >),
    (hpx_odeint_actions::operation_3_action< S , Operation >))

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename S , typename Operation >),
    (hpx_odeint_actions::operation_4_action< S , Operation >))

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename S >),
    (hpx_odeint_actions::identity_sync_action< S >))



struct dataflow_array_algebra
{

    //template< class S1 , class S2 , class S3 , class Op >
    // for now just a single state  type
    template< typename S , typename Op >
    void for_each3( S &s1 , const S &s2 , const S &s3 , Op op )
    {
        typedef hpx_odeint_actions::operation_3_action< typename S::result_type,Op> action;
        s1 = dataflow< action >( find_here() ,
                                 s1 , s2 , s3 , 
                                 op );
    }

    template< typename S , typename F1 , typename F2 >
    void for_each3( S &s1 , S &s2 , const S &s3 , dataflow_operations::scale_sum_swap2<F1,F2> op )
    {
        typedef hpx_odeint_actions::operation_3_action< typename S::result_type , dataflow_operations::scale_sum2<F1,F2> > action;
        typedef hpx_odeint_actions::identity_sync_action<typename S::result_type> copy_action;
        dataflow_base<typename S::result_type> tmp = s2;
        s2 = dataflow< action >( find_here() ,
                                 s1 , s2 , s3 , 
                                 dataflow_operations::scale_sum2<F1,F2>(op.m_alpha1,op.m_alpha2) );
        s1 = dataflow< copy_action >( find_here() , tmp , s2 );
    }

    template< typename S , typename Op >
    void for_each4( S &s1 , const S &s2 , const S &s3 , const S &s4 , Op op )
    {
        typedef hpx_odeint_actions::operation_4_action< typename S::result_type,Op> action;
        s1 = dataflow< action >( find_here() ,
                                 s1 , s2 , s3 , s4 ,
                                 op );
    }
};

#endif
