// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_OPERATIONS_HPP
#define DATAFLOW_SHARED_OPERATIONS_HPP

#include <vector>
#include <memory>

#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/runtime/actions/component_action.hpp>

struct dataflow_operations
{
    template< typename Fac1 , typename Fac2=Fac1 >
    struct scale_sum2
    {
        Fac1 m_alpha1;
        Fac2 m_alpha2;

        scale_sum2()
            : m_alpha1( 0 ) , m_alpha2( 0 )
        {}

        scale_sum2( Fac1 alpha1 , Fac2 alpha2 ) 
            : //managed_component_base() , 
              m_alpha1( alpha1 ) , m_alpha2( alpha2 ) 
        { }

        template< typename S1 , typename S2 , typename S3 >
        void operator() ( S1 &x1 , S2 &x2 , const S3 &x3 ) const
        {
            x1 = m_alpha1*x2 + m_alpha2*x3;
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // empty for now
        }
    };

    template< typename Fac1 , typename Fac2=Fac1 >
    struct scale_sum_swap2
    {
        Fac1 m_alpha1;
        Fac2 m_alpha2;

        scale_sum_swap2()
            : m_alpha1( 0 ) , m_alpha2( 0 )
        {}

        scale_sum_swap2( Fac1 alpha1 , Fac2 alpha2 ) 
            : //managed_component_base() , 
              m_alpha1( alpha1 ) , m_alpha2( alpha2 ) 
        { }

        template< typename S1 , typename S2 , typename S3 >
        void operator() ( S1 &x1 , S2 &x2 , const S3 &x3 ) const
        {
            const S1 tmp( x1 );
            x1 = m_alpha1 * x2 + m_alpha2 * x3;
            x2 = tmp;
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // empty for now
        }

    };


    template< typename Fac1 , typename Fac2=Fac1 , typename Fac3=Fac2>
    struct scale_sum3
    {
        Fac1 m_alpha1;
        Fac2 m_alpha2;
        Fac3 m_alpha3;

        scale_sum3()
            : m_alpha1( 0 ) , m_alpha2( 0 ) , m_alpha3( 0 )
        {}

        scale_sum3( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 ) 
            : //managed_component_base() , 
            m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) 
        { }

        template< typename S1 , typename S2 , typename S3 , typename S4 >
        void operator() ( S1 &x1 , const S2 &x2 , const S3 &x3 , const S4 &x4 ) const
        {
            x1 = m_alpha1*x2 + m_alpha2*x3 + m_alpha3*x4;
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // empty for now
        }
    };

    template< typename Fac1 , typename Fac2=Fac1 , typename Fac3=Fac2 , typename Fac4=Fac3 >
    struct scale_sum4
    {
        Fac1 m_alpha1;
        Fac2 m_alpha2;
        Fac3 m_alpha3;
        Fac4 m_alpha4;

        scale_sum4()
            : m_alpha1( 0 ) , m_alpha2( 0 ) , m_alpha3( 0 ) , m_alpha4( 0 )
        {}

        scale_sum4( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ) 
            : //managed_component_base() , 
            m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , 
            m_alpha3( alpha3 ) , m_alpha4( alpha4 )
        { }

        template< typename S1 , typename S2 , typename S3 , typename S4 , typename S5 >
        void operator() ( S1 &x1 , const S2 &x2 , const S3 &x3 , 
                          const S4 &x4 , const S5 &x5 ) const
        {
            x1 = m_alpha1*x2 + m_alpha2*x3 + m_alpha3*x4 + m_alpha4*x5;
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // empty for now
        }
    };

    template< typename Fac1 , typename Fac2=Fac1 , typename Fac3=Fac2 , 
              typename Fac4=Fac3 , typename Fac5=Fac4 >
    struct scale_sum5
    {
        Fac1 m_alpha1;
        Fac2 m_alpha2;
        Fac3 m_alpha3;
        Fac4 m_alpha4;
        Fac5 m_alpha5;

        scale_sum5()
            : m_alpha1( 0 ) , m_alpha2( 0 ) , m_alpha3( 0 ) , 
              m_alpha4( 0 ) , m_alpha5( 0 )
        {}

        scale_sum5( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 , Fac5 alpha5 ) 
            : //managed_component_base() , 
            m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , 
            m_alpha4( alpha4 ) , m_alpha5( alpha5 ) 
        { }

        template< typename S1 , typename S2 , typename S3 , typename S4 , 
                  typename S5 , typename S6 >
        void operator() ( S1 &x1 , const S2 &x2 , const S3 &x3 , 
                          const S4 &x4 , const S5 &x5 , const S6 &x6 ) const
        {
            x1 = m_alpha1*x2 + m_alpha2*x3 + m_alpha3*x4 + m_alpha4*x5 + 
                m_alpha5*x6;
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // empty for now
        }

    };


};

#endif
