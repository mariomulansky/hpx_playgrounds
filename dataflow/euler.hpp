// Copyright 2013 Mario Mulansky

#ifndef HPX_EULER_HPP
#define HPX_EULER_HPP

#include <vector>
#include <iostream>

#include <hpx/hpx.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>

//namespace hpx_step {

double euler_operation( double x , double dxdt , const double dt )
    {
        //std::clog << "euler operation" << std::endl;
        return x + dt*dxdt;
    }

//    HPX_DEFINE_PLAIN_ACTION(euler_operation, euler_operation_action);
//}

//HPX_REGISTER_PLAIN_ACTION( hpx_step::euler_operation_action );

HPX_PLAIN_ACTION( euler_operation, euler_operation_action );

namespace hpx_step {

    using hpx::lcos::dataflow;
    using hpx::lcos::dataflow_base;

    typedef dataflow_base<double> df_base;
    typedef std::vector< df_base > state_type;

    void euler_step( state_type &x , state_type &dxdt , const double dt , hpx::naming::id_type loc )
    {
        for( size_t i=0; i<boost::size(x) ; ++i )
        {
            //std::clog << "euler " << i << std::endl;
            x[i] = dataflow< euler_operation_action >( loc , x[i] , dxdt[i] , dt );
        }
    }
}

#endif
