#include <iostream>

//#define HPX_LIMIT 6

#include <hpx/hpx_main.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/iostreams.hpp>

#include "op.hpp"

using hpx::components::stub_base;
using hpx::components::client_base;
using hpx::components::managed_component;
using hpx::components::managed_component_base;
using hpx::components::simple_component;
using hpx::components::simple_component_base;

using hpx::find_here;
using hpx::async;
using hpx::lcos::future;

typedef boost_numeric_odeint::operations::my_op my_op;
typedef my_op::op_apply_action<double> op_apply_action;

int main()
{

    hpx::id_type target_locality = hpx::find_here();
    future<hpx::id_type> o = hpx::components::new_< my_op >( target_locality , 0.5 );

    op_apply_action act;

    act( o.get() , 2.0 );

    return 0;
}
