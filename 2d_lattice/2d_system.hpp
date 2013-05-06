// Copyright 2013 Mario Mulansky

#ifndef SYSTEM_2D_HPP
#define SYSTEM_2D_HPP

#include <vector>
#include <memory>
#include <cmath>

#include <boost/math/special_functions/sign.hpp>

#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/include/iostreams.hpp>

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

const double KAPPA = 3.5;
const double LAMBDA = 4.5;

namespace checked_math {
    inline double pow( double x , double y )
    {
        if( x==0.0 )
            // 0**y = 0, don't care for y = 0 or NaN
            return 0.0;
        using std::pow;
        return pow( x , y );
    }
}

double signed_pow( double x , double k )
{
    using boost::math::sign;
    using std::abs;
    return checked_math::pow( abs(x) , k ) * sign(x);
}

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vecvec;
typedef std::vector< dataflow_base< shared_vec > > state_type;


shared_vecvec system_first_block( shared_vecvec q , shared_vecvec q_d , shared_vecvec dpdt )
{
    const size_t N = q->size();

    //hpx::cout << (boost::format("system first block sizes %d , %d\n") % N % M) << hpx::flush;

    double coupling_lr = 0.0;
    const size_t M = (*q)[0].size();
    dvec coupling_ud( M , 0.0 );
    for( size_t i=0 ; i<N-1 ; ++i )
    {
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[i][j] = -signed_pow( (*q)[i][j] , KAPPA-1 ) 
                + coupling_lr + coupling_ud[j];
            coupling_lr = signed_pow( (*q)[i][j]-(*q)[i][j+1] , LAMBDA-1 );
            coupling_ud[j] = signed_pow( (*q)[i][j]-(*q)[i+1][j] , LAMBDA-1 );
            (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[i][M-1] = -signed_pow( (*q)[i][M-1] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[M-1];
        coupling_ud[M-1] = signed_pow( (*q)[i][M-1]-(*q)[i+1][M-1] , LAMBDA-1 );
        coupling_lr = 0.0;
        (*dpdt)[i][M-1] -= coupling_ud[M-1];
    }
    for( size_t j=0 ; j<M-1 ; ++j )
    {
        (*dpdt)[N-1][j] = -signed_pow( (*q)[N-1][j] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[j];
        coupling_lr = signed_pow( (*q)[N-1][j]-(*q)[N-1][j+1] , LAMBDA-1 );
        (*dpdt)[N-1][j] -= coupling_lr + signed_pow( (*q)[N-1][j] - (*q_d)[0][j] , LAMBDA-1 );
    }
    (*dpdt)[N-1][M-1] = -signed_pow( (*q)[N-1][M-1] , KAPPA-1 ) 
        + coupling_lr + coupling_ud[M-1]
        - signed_pow( (*q)[N-1][M-1] - (*q_d)[0][M-1] , LAMBDA-1 );
    
    return dpdt;
}

HPX_PLAIN_ACTION(system_first_block, system_first_block_action);

shared_vecvec system_center_block( shared_vecvec q , shared_vecvec q_u , 
                                   shared_vecvec q_d , shared_vecvec dpdt )
{
    using checked_math::pow;
    const size_t N = q->size();
    const size_t M = (*q)[0].size();

    //hpx::cout << (boost::format("system center block sizes %d , %d\n") % N % M) << hpx::flush;

    double coupling_lr = 0.0;
    dvec coupling_ud( M );

    for( size_t j=0 ; j<M-1 ; ++j )
    {
        (*dpdt)[0][j] = -signed_pow( (*q)[0][j] , KAPPA-1 ) + coupling_lr
            - signed_pow( (*q)[0][j] - (*q_u)[q_u->size()-1][j] , LAMBDA-1 );
        coupling_lr = signed_pow( (*q)[0][j]-(*q)[0][j+1] , LAMBDA-1 );
        coupling_ud[j] = signed_pow( (*q)[0][j]-(*q)[1][j] , LAMBDA-1 );
        (*dpdt)[0][j] -= coupling_lr + coupling_ud[j];
    }
    (*dpdt)[0][M-1] = -signed_pow( (*q)[0][M-1] , KAPPA-1 ) + coupling_lr
        - signed_pow( (*q)[0][M-1]-(*q_u)[q_u->size()-1][M-1] , LAMBDA-1 );
    coupling_ud[M-1] = signed_pow( (*q)[0][M-1]-(*q)[1][M-1] , LAMBDA-1 );
    coupling_lr = 0.0;
    (*dpdt)[0][M-1] -= coupling_ud[M-1];


    for( size_t i=1 ; i<N-1 ; ++i )
    {
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[i][j] = -signed_pow( (*q)[i][j] , KAPPA-1 ) 
                + coupling_lr + coupling_ud[j];
            coupling_lr = signed_pow( (*q)[i][j]-(*q)[i][j+1] , LAMBDA-1 );
            coupling_ud[j] = signed_pow( (*q)[i][j]-(*q)[i+1][j] , LAMBDA-1 );
            (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[i][M-1] = -signed_pow( (*q)[i][M-1] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[M-1];
        coupling_ud[M-1] = signed_pow( (*q)[i][M-1]-(*q)[i+1][M-1] , LAMBDA-1 );
        coupling_lr = 0.0;
        (*dpdt)[i][M-1] -= coupling_ud[M-1];
    }

    for( size_t j=0 ; j<M-1 ; ++j )
    {
        (*dpdt)[N-1][j] = -signed_pow( (*q)[N-1][j] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[j];
        coupling_lr = signed_pow( (*q)[N-1][j]-(*q)[N-1][j+1] , LAMBDA-1 );
        (*dpdt)[N-1][j] -= coupling_lr + signed_pow( (*q)[N-1][j] - (*q_d)[0][j] , LAMBDA-1 );
    }
    (*dpdt)[N-1][M-1] = -signed_pow( (*q)[N-1][M-1] , KAPPA-1 ) 
        + coupling_lr + coupling_ud[M-1]
        - signed_pow( (*q)[N-1][M-1] - (*q_d)[0][M-1] , LAMBDA-1 );
    
    return dpdt;
}

HPX_PLAIN_ACTION(system_center_block, system_center_block_action);


shared_vecvec system_last_block( shared_vecvec q , shared_vecvec q_u , shared_vecvec dpdt )
{
    using checked_math::pow;
    const size_t N = q->size();
    const size_t M = (*q)[0].size();

    //hpx::cout << (boost::format("system last block sizes %d , %d\n") % N % M) << hpx::flush;

    double coupling_lr = 0.0;
    dvec coupling_ud( M );

    for( size_t j=0 ; j<M-1 ; ++j )
    {
        (*dpdt)[0][j] = -signed_pow( (*q)[0][j] , KAPPA-1 ) + coupling_lr
            - signed_pow( (*q)[0][j] - (*q_u)[q_u->size()-1][j] , LAMBDA-1 );
        coupling_lr = signed_pow( (*q)[0][j]-(*q)[0][j+1] , LAMBDA-1 );
        coupling_ud[j] = signed_pow( (*q)[0][j]-(*q)[1][j] , LAMBDA-1 );
        (*dpdt)[0][j] -= coupling_lr + coupling_ud[j];
    }
    (*dpdt)[0][M-1] = -signed_pow( (*q)[0][M-1] , KAPPA-1 ) + coupling_lr
        - signed_pow( (*q)[0][M-1]-(*q_u)[q_u->size()-1][M-1] , LAMBDA-1 );
    coupling_ud[M-1] = signed_pow( (*q)[0][M-1]-(*q)[1][M-1] , LAMBDA-1 );
    coupling_lr = 0.0;
    (*dpdt)[0][M-1] -= coupling_ud[M-1];


    for( size_t i=1 ; i<N-1 ; ++i )
    {
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[i][j] = -signed_pow( (*q)[i][j] , KAPPA-1 ) 
                + coupling_lr + coupling_ud[j];
            coupling_lr = signed_pow( (*q)[i][j]-(*q)[i][j+1] , LAMBDA-1 );
            coupling_ud[j] = signed_pow( (*q)[i][j]-(*q)[i+1][j] , LAMBDA-1 );
            (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[i][M-1] = -signed_pow( (*q)[i][M-1] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[M-1];
        coupling_ud[M-1] = signed_pow( (*q)[i][M-1]-(*q)[i+1][M-1] , LAMBDA-1 );
        coupling_lr = 0.0;
        (*dpdt)[i][M-1] -= coupling_ud[M-1];
    }

    for( size_t j=0 ; j<M-1 ; ++j )
    {
        (*dpdt)[N-1][j] = -signed_pow( (*q)[N-1][j] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[j];
        coupling_lr = signed_pow( (*q)[N-1][j]-(*q)[N-1][j+1] , LAMBDA-1 );
        (*dpdt)[N-1][j] -= coupling_lr;
    }
    (*dpdt)[N-1][M-1] = -signed_pow( (*q)[N-1][M-1] , KAPPA-1 ) 
        + coupling_lr + coupling_ud[M-1];
    
    return dpdt;
}

HPX_PLAIN_ACTION(system_last_block, system_last_block_action);

void system_2d( const state_type &q , state_type &dpdt )
{
    const size_t N = q.size();
    // first row
    dpdt[0] = dataflow< system_first_block_action >( find_here() , q[0] , q[1] , dpdt[0] );
    // middle rows
    for( size_t i=1 ; i<N-1 ; i++ )
    {
        dpdt[i] = dataflow< system_center_block_action >( find_here() , q[i] , q[i-1] , q[i+1] , dpdt[i] );
    }
    dpdt[N-1] = dataflow< system_last_block_action >( find_here() , q[N-1] , q[N-2] , dpdt[N-1] );
}


double energy( const dvecvec &q , const dvecvec &p )
{
    using checked_math::pow;
    using std::abs;
    const size_t N = q.size();
    double energy = 0.0;
    for( size_t i=0 ; i<N-1 ; ++i )
    {
        const size_t M = q[i].size();
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            energy += 0.5*p[i][j]*p[i][j] + pow( q[i][j] , KAPPA ) / KAPPA
                + pow( abs(q[i][j]-q[i][j+1]) , LAMBDA ) / LAMBDA
                + pow( abs(q[i][j]-q[i+1][j]) , LAMBDA ) / LAMBDA;
        }
        energy += 0.5*p[i][M-1]*p[i][M-1] + pow( q[i][M-1] , KAPPA ) / KAPPA
            + pow( abs(q[i][M-1]-q[i+1][M-1]) , LAMBDA ) / LAMBDA;
    }
    const size_t M = q[N-1].size();
    for( size_t j=0 ; j<M-1 ; ++j )
    {
        energy += 0.5*p[N-1][j]*p[N-1][j] + pow( q[N-1][j] , KAPPA ) / KAPPA
            + pow( abs(q[N-1][j]-q[N-1][j+1]) , LAMBDA ) / LAMBDA;
    }
    energy += 0.5*p[N-1][M-1]*p[N-1][M-1] + pow( q[N-1][M-1] , KAPPA ) / KAPPA;
    return energy;
}

template< typename S >
double energy( const S &q_fut , const S &p_fut )
{
    dvecvec q,p;
    for( size_t i=0 ; i<q_fut.size() ; ++i )
    {
        for( size_t j=0 ; j<(*(q_fut[i].get())).size() ; ++j )
        {
            q.push_back( (*(q_fut[i].get()))[j] );
            p.push_back( (*(p_fut[i].get()))[j] );
        }
    }
    return energy( q , p );
}

#endif
