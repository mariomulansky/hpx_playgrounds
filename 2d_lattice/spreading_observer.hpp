// Copyright Mario Mulansky 2013
#ifndef SPREADING_OBSERVER_HPP
#define SPREADING_OBSERVER_HPP

#include <vector>
#include <utility>

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/include/iostreams.hpp>

#include "hpx_odeint_actions.hpp"
#include "2d_system.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

typedef dataflow_base< double > df_base;
typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;
typedef std::vector< dataflow_base< shared_vec > > state_type;

const double border = 1E-50;

double excitation_area( shared_vec q , shared_vec p )
{
    hpx::cout << "excitation area\n" << hpx::flush;

    const int N = q->size();
    const int M = (*q)[0].size();

    double a = 0.0;

    for( int i=0 ; i<N ; ++i )
        for( int j=0 ; j<M ; ++j )
        {
            double loc_e = 0.5*(*p)[i][j]*(*p)[i][j] + checked_math::pow( (*q)[i][j] , KAPPA )/KAPPA;
            if( loc_e > border )
                a += 1.0;
        }

    return a;
}

HPX_PLAIN_ACTION( excitation_area , excitation_area_action );

struct spreading_observer
{

    std::vector< std::pair< double , std::vector< dataflow_base<double> > > > m_values;

    void operator()( state_type &q , state_type &p , double t )
    {
        int N = q.size();
        std::vector< dataflow_base<double> > v(N);
        for( int i=0 ; i<N ; ++i )
        {
            v[i] = dataflow< excitation_area_action >( find_here() ,
                                                       q[i] ,
                                                       p[i] );
            // synchronize
            p[i] = dataflow< sync1_action<shared_vec,double> >( find_here() , p[i] , v[i] );
            q[i] = dataflow< sync1_action<shared_vec,double> >( find_here() , q[i] , v[i] );
        }
        m_values.push_back( make_pair( t , v ) );
    }

};

#endif
