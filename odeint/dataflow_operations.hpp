// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_OPERATIONS_HPP
#define DATAFLOW_OPERATIONS_HPP

#include <vector>

#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/runtime/actions/component_action.hpp>

typedef std::vector<double> dvec;

using hpx::components::simple_component;
using hpx::components::simple_component_base;
using hpx::components::managed_component;
using hpx::components::managed_component_base;

// double hpx_scale_sum2( double x , double dxdt , double a1 , double a2 )
// {
//     //std::clog << "euler operation" << std::endl;
//     return a1*x + a2*dxdt;
// }

// HPX_PLAIN_ACTION( hpx_scale_sum2, hpx_scale_sum2_action );

// dvec hpx_scale_sum2_sub( dvec out , const dvec &x , const dvec &dxdt , double a1 , double a2 )
// {
//     //std::clog << "euler operation" << std::endl;
//     for( size_t i=0 ; i<x.size() ; ++i )
//         out[i] = a1*x[i] + a2*dxdt[i];
//     return out;
// }

// HPX_PLAIN_ACTION( hpx_scale_sum2_sub, hpx_scale_sum2_sub_action );

/*dvec_ref hpx_scale_sum2_ref( devc_ref x_ , dvec_ref dxdt_ , double a1 , double a2 )
{
    dvec &x = x_;
    dvec &dxdt = dxdt_;
    for( int i=0 ; i<x.size() ; ++i )
        x[i] += a2*dxdt[i];
    return x_;
}

HPX_PLAIN_ACTION( hpx_scale_sum2_ref, hpx_scale_sum2_ptr_action );
*/


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


struct dataflow_actions
{
    //typedef dvec vector_type;

    struct operation3 : managed_component_base< operation3 >
    {

        template< typename Vector , typename Operation >
        Vector apply( Vector x1 , const Vector &x2 , const Vector &x3 , const Operation &op )
        {
            for( size_t i=0 ; i<x1.size() ; ++i )
                op( x1[i] , x2[i] , x3[i] );
            return x1;
        }

        // action definition for template member function
        template <typename Vector , typename Operation>
        struct operation3_action
            : hpx::actions::make_action<Vector (operation3::*)( Vector , const Vector & , const Vector & , const Operation& ),
                                        &operation3::template apply< Vector , Operation >, operation3_action< Vector , Operation > >
        {};
        
    };

};

HPX_REGISTER_ACTION_DECLARATION_TEMPLATE(
                                         (template <typename Vector , typename Operation>),
                                         (dataflow_actions::operation3::operation3_action<Vector , Operation>)
)


#endif
