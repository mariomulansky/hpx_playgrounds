#include <hpx/hpx.hpp>
#include <hpx/runtime/components/component_factory.hpp>

#include <hpx/util/portable_binary_iarchive.hpp>
#include <hpx/util/portable_binary_oarchive.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/export.hpp>

#include "dataflow_operations.hpp"

using hpx::components::managed_component;

typedef dataflow_actions::operation3 operation3;
typedef managed_component< operation3 > operation3_component;

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY( operation3_component , operation3 );

// important!
// typedef scale_sum2::scale_sum2_action scale_sum2_action;

// HPX_REGISTER_ACTION( scale_sum2_action );
