// Copyright 2012 Mario Mulansky

#ifndef SPREADING_OBSERVER_HPP
#define SPREADING_OBSERVER_HPP

#include <vector>
#include <utility>

#include <boost/numeric/odeint/util/unwrap_reference.hpp>

struct spreading_observer
{


    const double m_beta;
    const double m_kap;
    const double m_lam;
    const int m_block_size;
    int m_t;
    std::vector< std::pair< double , double > > m_values;

    spreading_observer( const double kap , const double lam , 
                        const double beta , const int block_size = 1 )
        : m_kap( kap ) , m_lam( lam ) , m_beta( beta ) , 
          m_block_size( block_size ) , m_t(0)
    { }


    template< typename State >
    void operator()( const State &x , double t )
    {
        if( (m_t%10) == 0 )
        {
            using boost::numeric::odeint::unwrap_reference;
            typedef typename unwrap_reference< State >::type state_type;
            typedef typename unwrap_reference< typename state_type::first_type >::type coor_in_type;
            typedef typename unwrap_reference< typename state_type::second_type >::type momentum_in_type;
            const state_type &state = x;
            const coor_in_type &q = state.first;
            const momentum_in_type &p = state.second;
        
            using checked_math::pow;
            // q and dpdt are 2d
            const size_t N = q.size();
            const size_t M = q[0].size();

            double a = 0.0;

#ifndef NO_OMP
#pragma omp parallel for reduction(+:a) schedule( runtime )
#endif //NO_OMP
            for( size_t i=0 ; i<N-1 ; ++i )
            {
                for( size_t j=0 ; j<M-1 ; ++j )
                {
                    double loc_en = p[i][j]*p[i][j] / 2.0 + checked_math::pow( q[i][j] , m_kap ) / m_kap;
                    if( loc_en > 1E-50 )
                    {
                        a += 1;
                    }
                }
            }
            m_values.push_back( std::make_pair( t , a ) );
        }
        m_t++;
    }

};

#endif
