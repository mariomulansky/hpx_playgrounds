#include <boost/mpl/print.hpp>

#include <hpx/hpx.hpp>
#include <hpx/runtime/components/component_factory.hpp>

#include <hpx/util/portable_binary_iarchive.hpp>
#include <hpx/util/portable_binary_oarchive.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/export.hpp>

#include "op.hpp"

using hpx::components::simple_component;
using hpx::components::simple_component_base;

typedef boost_numeric_odeint::operations::my_op my_op;
typedef simple_component< my_op > op_type;

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY( op_type, my_op );

// typedef my_op::op_apply_action op_apply_action;

// HPX_REGISTER_ACTION( op_apply_action );
