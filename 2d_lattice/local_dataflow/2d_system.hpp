// Copyright 2013 Mario Mulansky

#ifndef SYSTEM_2D_HPP
#define SYSTEM_2D_HPP

#include <vector>
#include <memory>
#include <cmath>

#include <boost/math/special_functions/sign.hpp>
#include <boost/thread/thread.hpp>

#include <hpx/lcos/future.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/util/unwrap.hpp>

using hpx::lcos::local::dataflow;
using hpx::lcos::future;
using hpx::lcos::wait;
using hpx::util::unwrap;

const double KAPPA = 3.5;
const double LAMBDA = 4.5;

namespace checked_math {
    inline double pow( double x , double y )
    {
        if( x==0.0 )
            // 0**y = 0, don't care for y = 0 or NaN
            return 0.0;
        using std::pow;
        using std::abs;
        return pow( abs(x) , y );
    }
}

double signed_pow( double x , double k )
{
    using boost::math::sign;
    using std::abs;
    return checked_math::pow( x , k ) * sign(x);
}

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vecvec;
typedef std::vector< future< shared_vec > > state_type;

struct system_first_block
{
    shared_vecvec operator()( shared_vecvec q , const dvec q_d , shared_vecvec dpdt ) const
    {
        //hpx::cout << (boost::format("first block\n") ) << hpx::flush;

        const size_t N = q->size();

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
            (*dpdt)[N-1][j] -= coupling_lr + signed_pow( (*q)[N-1][j] - q_d[j] , LAMBDA-1 );
        }
        (*dpdt)[N-1][M-1] = -signed_pow( (*q)[N-1][M-1] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[M-1]
            - signed_pow( (*q)[N-1][M-1] - q_d[M-1] , LAMBDA-1 );
    
        //hpx::cout << (boost::format("first block done\n") ) << hpx::flush;
        return dpdt;
    }
};


struct system_center_block
{

    shared_vecvec operator() ( shared_vecvec q , const dvec q_u , 
                               const dvec q_d , shared_vecvec dpdt ) const
    {
        //hpx::cout << (boost::format("center block\n") ) << hpx::flush;
 
        using checked_math::pow;
        const size_t N = q->size();
        const size_t M = (*q)[0].size();

        double coupling_lr = 0.0;
        dvec coupling_ud( M );

        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[0][j] = -signed_pow( (*q)[0][j] , KAPPA-1 ) + coupling_lr
                - signed_pow( (*q)[0][j] - q_u[j] , LAMBDA-1 );
            coupling_lr = signed_pow( (*q)[0][j]-(*q)[0][j+1] , LAMBDA-1 );
            coupling_ud[j] = signed_pow( (*q)[0][j]-(*q)[1][j] , LAMBDA-1 );
            (*dpdt)[0][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[0][M-1] = -signed_pow( (*q)[0][M-1] , KAPPA-1 ) + coupling_lr
            - signed_pow( (*q)[0][M-1]-q_u[M-1] , LAMBDA-1 );
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
            (*dpdt)[N-1][j] -= coupling_lr + signed_pow( (*q)[N-1][j] - q_d[j] , LAMBDA-1 );
        }
        (*dpdt)[N-1][M-1] = -signed_pow( (*q)[N-1][M-1] , KAPPA-1 ) 
            + coupling_lr + coupling_ud[M-1]
            - signed_pow( (*q)[N-1][M-1] - q_d[M-1] , LAMBDA-1 );
    
        //hpx::cout << (boost::format("center block done\n") ) << hpx::flush;

        return dpdt;
    }
};


struct system_last_block
{

    typedef shared_vec result_type;

    shared_vecvec operator()( shared_vecvec q , const dvec q_u , shared_vecvec dpdt ) const
    {
        //hpx::cout << (boost::format("last block\n") ) << hpx::flush;

        using checked_math::pow;
        const size_t N = q->size();
        const size_t M = (*q)[0].size();

        double coupling_lr = 0.0;
        dvec coupling_ud( M );

        //hpx::cout << (boost::format("last block iterating...\n") ) << hpx::flush;

        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[0][j] = -signed_pow( (*q)[0][j] , KAPPA-1 ) + coupling_lr
                - signed_pow( (*q)[0][j] - q_u[j] , LAMBDA-1 );
            coupling_lr = signed_pow( (*q)[0][j]-(*q)[0][j+1] , LAMBDA-1 );
            coupling_ud[j] = signed_pow( (*q)[0][j]-(*q)[1][j] , LAMBDA-1 );
            (*dpdt)[0][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[0][M-1] = -signed_pow( (*q)[0][M-1] , KAPPA-1 ) + coupling_lr
            - signed_pow( (*q)[0][M-1]-q_u[M-1] , LAMBDA-1 );
        coupling_ud[M-1] = signed_pow( (*q)[0][M-1]-(*q)[1][M-1] , LAMBDA-1 );
        coupling_lr = 0.0;
        (*dpdt)[0][M-1] -= coupling_ud[M-1];

        for( size_t i=1 ; i<N-1 ; ++i )
        {
            //hpx::cout << (boost::format("row %d\n") % i ) << hpx::flush;
            for( size_t j=0 ; j<M-1 ; ++j )
            {
                //hpx::cout << (boost::format("row %d , col %d \n") % i %j ) << hpx::flush;
                (*dpdt)[i][j] = -signed_pow( (*q)[i][j] , KAPPA-1 ) 
                    + coupling_lr + coupling_ud[j];
                coupling_lr = signed_pow( (*q)[i][j]-(*q)[i][j+1] , LAMBDA-1 );
                coupling_ud[j] = signed_pow( (*q)[i][j]-(*q)[i+1][j] , LAMBDA-1 );
                (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
            }
            //hpx::cout << (boost::format("row %d , col %d \n") % i % (M-1) ) << hpx::flush;
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
    
        //hpx::cout << (boost::format("last block done\n") ) << hpx::flush;

        return dpdt;
    }
};

void system_2d( state_type &q , state_type &dpdt )
{
    // works on shared data, but coupling data is provided as copy
    const size_t N = q.size();

    //hpx::cout << boost::format("system call size: %d , %d ...\n") % (q.size()) % (dpdt.size()) << hpx::flush;

    //state_type dpdt_(N);
    // first row
    dpdt[0] = dataflow( hpx::launch::async , unwrap(system_first_block()) , q[0] , 
                        dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
        { return (*v)[0]; }) , q[1] ) , 
                        dpdt[0] );
    // middle rows
    for( size_t i=1 ; i<N-1 ; i++ )
    {
        dpdt[i] = dataflow( hpx::launch::async , unwrap(system_center_block()) , q[i] , 
                            dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
            { return (*v)[v->size()-1]; }) ,
                                      q[i-1] ) , 
                            dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
            { return (*v)[0]; }) , q[i+1] ) ,
                            dpdt[i] );
    }
    dpdt[N-1] = dataflow( hpx::launch::async , unwrap(system_last_block()) , q[N-1] , 
                          dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
        { return (*v)[v->size()-1]; }), 
                                    q[N-2] ) , 
                          dpdt[N-1] );
    
    //hpx::cout << "system call finished\n" << hpx::flush;

    /*
    // coupling synchronization step
    dpdt[0] = dataflow( hpx::launch::sync , 
                        unwrap([]( shared_vec x , shared_vec sync ){ return x; } ) , 
                        dpdt[0] , dpdt[1] );

    for( size_t i=1 ; i<N-1 ; i++ )
    {
        dpdt[i] = dataflow( hpx::launch::sync , 
                            unwrap([]( shared_vec x , shared_vec sync1 , shared_vec sync2 )
            { return x; } ) , 
                            dpdt[i] , dpdt[i-1] , dpdt[i+1] );
    }
    dpdt[N-1] = dataflow( hpx::launch::sync , 
                           unwrap([]( shared_vec x , shared_vec sync ){ return x; } ) , 
                           dpdt[N-1] , dpdt[N-2] );
    */
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
