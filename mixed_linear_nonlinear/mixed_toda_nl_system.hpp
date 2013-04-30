// Copyright 2013 Mario Mulansky
//
#ifndef MIXED_TODA_NL_SYSTEM_HPP
#define MIXED_TODA_NL_SYSTEM_HPP

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/sign.hpp>

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::lcos::future;

typedef dataflow_base< double > df_base;
typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef std::vector< dataflow_base< shared_vec > > state_type;
typedef std::vector< dvec > vec_dvec;

const int KAPPA=3.5;
const int LAMBDA=4.5;


// left part of the chain (Toda)
shared_vec osc_op1( shared_vec q , shared_vec q_r , shared_vec dpdt )
{
    using std::pow;
    using std::abs;
    using std::exp;
    using boost::math::sign;

    const size_t M = q->size();

    (*dpdt)[0] = - exp( (*q)[0] - (*q)[1] );
    for( size_t i=1 ; i<M-1 ; ++i )
    {
        (*dpdt)[i] = exp( (*q)[i-1] - (*q)[i] ) - exp( (*q)[i] - (*q)[i+1] );
    }
    (*dpdt)[M-1] = exp( (*q)[M-2] - (*q)[M-1] )
        // power-law coupling to the center part
        - pow( abs((*q)[M-1] - (*q_r)[0]) , LAMBDA-1 ) * sign((*q)[M-1]-(*q_r)[0]);

    return dpdt;
}

HPX_PLAIN_ACTION(osc_op1, osc_op1_action);

// center part of the chain (power laws)
shared_vec osc_op2( shared_vec q , shared_vec q_l , shared_vec q_r , shared_vec dpdt )
{
    using std::pow;
    using std::abs;
    using boost::math::sign;

    const size_t M = q->size();

    double coupling = pow( abs((*q_l)[q_l->size()-1]-(*q)[0]) , LAMBDA-1 ) * sign((*q_l)[q_l->size()-1]-(*q)[0]);

    for( size_t i=0 ; i<M-1 ; ++i )
    {
        (*dpdt)[i] = -pow( abs((*q)[i]) , KAPPA-1 ) * sign((*q)[i])
            + coupling;
        coupling = pow( abs((*q)[i]-(*q)[i+1]) , LAMBDA-1 ) * sign((*q)[i]-(*q)[i+1]);
        (*dpdt)[i] -= coupling;
    }

    (*dpdt)[M-1] = -pow( abs((*q)[M-1]) , KAPPA-1 ) * sign((*q)[M-1]) 
        - pow( abs((*q)[M-1] - (*q_r)[0]) , LAMBDA-1 ) * sign((*q)[M-1]-(*q_r)[0])
        + coupling;

    return dpdt;
}

HPX_PLAIN_ACTION(osc_op2, osc_op2_action);

// right part of the chain (linear)
shared_vec osc_op3( shared_vec q , shared_vec q_l , shared_vec dpdt )
{
    using std::pow;
    using std::abs;
    using boost::math::sign;

    const size_t M = q->size();

    // left coupling, nonlinear
    (*dpdt)[0] = -exp( (*q)[0] - (*q)[1] )
        + pow( abs((*q_l)[M-1] - (*q)[0]) , LAMBDA-1 ) * sign((*q_l)[M-1]-(*q)[0]);

    for( size_t i=1 ; i<M-1 ; ++i )
    {
        (*dpdt)[i] = exp( (*q)[i-1] - (*q)[i] ) - exp( (*q)[i] - (*q)[i+1] );
    }

    // no coupling to the right
    (*dpdt)[M-1] = exp( (*q)[M-2] - (*q)[M-1] );

    return dpdt;
}

HPX_PLAIN_ACTION(osc_op3, osc_op3_action);


void mixed_toda_nl_system( const state_type &q , state_type &dpdt )
{
    //std::cout << "system call with N=" << N << std::endl;
    dpdt[0] = dataflow< osc_op1_action >( find_here() ,
                                          q[0] , q[1] ,
                                          dpdt[0] );

    dpdt[1] = dataflow< osc_op2_action >( find_here() , 
                                          q[1] , q[0] , q[2] ,
                                          dpdt[1] );

    dpdt[2] = dataflow< osc_op3_action >( find_here() ,
                                          q[2] , q[1] ,
                                          dpdt[2] );
    //std::cout << "-------" << std::endl;
}

double energy( const dvec &q1 , const dvec &q2 , const dvec &q3 ,
               const dvec &p1 , const dvec &p2 , const dvec &p3 )
{
    using std::pow;
    using std::abs;
    using boost::math::pow;

    const size_t N1=q1.size();
    double e(0.0);
    // first part (linear)
    for( size_t n=0 ; n<N1-1 ; ++n )
    {
        e += 0.5*pow<2>( p1[n] ) + exp( q1[n] - q1[n+1] );
    }
    // nl coupling
    e += 0.5*pow<2>( p1[N1-1] ) + pow( abs(q1[N1-1]-q2[0]) , LAMBDA) / LAMBDA;

    const size_t N2=q2.size();
    // second part (nonlinear)
    for( size_t n=0 ; n<N2-1 ; ++n )
    {
        e += 0.5*pow<2>( p2[n] )
            + pow( abs(q2[n]) , KAPPA )/KAPPA 
            + pow( abs(q2[n]-q2[n+1]) , LAMBDA )/LAMBDA;
    }
    //nl coupling
    e += 0.5*pow<2>( p2[N2-1] ) 
        + pow( abs(q2[N2-1]) , KAPPA )/KAPPA
        + pow( abs(q2[N2-1]-q3[0]) , LAMBDA )/LAMBDA;

    const size_t N3=q3.size();
    // third part (linear)
    for( size_t n=0 ; n<N3-1 ; ++n )
    {
        e += 0.5*pow<2>( p3[n] ) + exp( q3[n] - q3[n+1] );
    }
    e += 0.5*pow<2>( p3[N3-1] );

    return e;
}

double energy( const std::vector< future<shared_vec> > &q_futures ,
               const std::vector< future<shared_vec> > &p_futures )
{
    dvec q1( q_futures[0].get()->begin() , q_futures[0].get()->end() );
    dvec q2( q_futures[1].get()->begin() , q_futures[1].get()->end() );
    dvec q3( q_futures[2].get()->begin() , q_futures[2].get()->end() );
    dvec p1( p_futures[0].get()->begin() , p_futures[0].get()->end() );
    dvec p2( p_futures[1].get()->begin() , p_futures[1].get()->end() );
    dvec p3( p_futures[2].get()->begin() , p_futures[2].get()->end() );

    return energy( q1 , q2 , q3 , p1 , p2 , p3 );
}

#endif
