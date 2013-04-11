// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_OPERATIONS_HPP
#define DATAFLOW_OPERATIONS_HPP

#include <vector>

typedef std::vector<double> dvec;

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

double hpx_scale_sum2( double x , double dxdt , double a1 , double a2 )
{
    //std::clog << "euler operation" << std::endl;
    return a1*x + a2*dxdt;
}

HPX_PLAIN_ACTION( hpx_scale_sum2, hpx_scale_sum2_action );

dvec hpx_scale_sum2_sub( dvec x , dvec dxdt , double a1 , double a2 )
{
    //std::clog << "euler operation" << std::endl;
    dvec out( x.size() );
    for( size_t i=0 ; i<x.size() ; ++i )
        out[i] = a1*x[i] + a2*dxdt[i];
    return out;
}

HPX_PLAIN_ACTION( hpx_scale_sum2_sub, hpx_scale_sum2_sub_action );

struct dataflow_operations
{

    template< class Fac1 = double , class Fac2 = Fac1 >
    struct scale_sum2
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum2( Fac1 alpha1 , Fac2 alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) 
        { }
    };

};


#endif
