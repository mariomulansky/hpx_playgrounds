#ifndef OP_HPP
#define OP_HPP

#include <iostream>

//#include <hpx/include/components.hpp>
//#include <hpx/include/actions.hpp>
#include <hpx/runtime/components/server/simple_component_base.hpp>
#include <hpx/runtime/actions/component_action.hpp>


using hpx::components::simple_component;
using hpx::components::simple_component_base;

namespace boost_numeric_odeint {

    struct operations {

        struct my_op : simple_component_base< my_op >
        {
            my_op() : m_a( 1.0 )
            { }


            my_op( double a ) : m_a( a )
            { }

            template< typename T >
            T apply( T x )
            {
                std::cout << "apply " << m_a << " , " << m_a*x << std::endl;
                return m_a*x;
            }

            template <typename T>
            struct op_apply_action
                : hpx::actions::make_action<T (my_op::*)(T),
                                            &my_op::template apply<T>, op_apply_action<T> >
            {};

            //HPX_DEFINE_COMPONENT_ACTION( my_op , apply , op_apply_action );
            
            double m_a;
        };

    };

}

HPX_REGISTER_ACTION_DECLARATION_TEMPLATE(
    (template <typename T>),
    (boost_numeric_odeint::operations::my_op::op_apply_action<T>)
)

// HPX_REGISTER_ACTION_DECLARATION( boost_numeric_odeint::operations::my_op::op_apply_action , 
//                                  op_apply_action );

#endif
