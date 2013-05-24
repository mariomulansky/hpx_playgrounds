// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_OPERATIONS_HPP
#define DATAFLOW_SHARED_OPERATIONS_HPP

#include <vector>
#include <memory>

#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/runtime/actions/component_action.hpp>

#include "serialization.hpp"

typedef std::vector<double> dvec;

using hpx::components::managed_component;
using hpx::components::managed_component_base;

struct dataflow_operations
{
    template< typename Fac1 , typename Fac2=Fac1 >
    struct scale_sum2
    {
        Fac1 m_alpha1;
        Fac2 m_alpha2;

        scale_sum2()
            : m_alpha1( 0 ) , m_alpha2( 0 )
        {}

        scale_sum2( Fac1 alpha1 , Fac2 alpha2 ) 
            : //managed_component_base() , 
              m_alpha1( alpha1 ) , m_alpha2( alpha2 ) 
        { }

        template< typename S1 , typename S2 , typename S3 >
        void operator() ( S1 &x1 , const S2 &x2 , const S3 &x3 ) const
        {
            x1 = m_alpha1*x2 + m_alpha2*x3;
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // empty for now
        }

    };

};

// namespace hpx_odeint_actions {

//     template< typename S , typename Operation >
//     S operation3( S x1 , const S &x2 , const S &x3 , 
//                   Operation op )
//     {
//         for( size_t i=0 ; i<x1->size() ; ++i )
//             op( (*x1)[i] , (*x2)[i] , (*x3)[i] );
//         return x1;
//     }


//     template< typename S , typename Operation >
//     struct operation3_action
//         : hpx::actions::make_action<
//         S (*)( S , 
//                const S& , 
//                const S& , 
//                Operation op ) , 
//         &operation3<S,Operation>, 
//         operation3_action<S,Operation> >
//     {};

// }

// HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
//     (template< typename S , typename Operation >),
//     (hpx_odeint_actions::operation3_action< S , Operation >))

#endif
