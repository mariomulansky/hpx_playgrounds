// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_OPERATIONS2D_HPP
#define DATAFLOW_SHARED_OPERATIONS2D_HPP

#include <vector>
#include <memory>

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;

struct dataflow_shared_operations2d
{
    template< typename Fac1 , typename Fac2=Fac1 >
    struct scale_sum2
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum2()
            : m_alpha1( 0 ) , m_alpha2( 0 )
        {}

        scale_sum2( Fac1 alpha1 , Fac2 alpha2 ) 
              m_alpha1( alpha1 ) , m_alpha2( alpha2 ) 
        { }

        template< typename S1 , typename S2 , typename S3 >
        void operator() ( S1 &x1 , const S2 &x2 , const S3 &x3 ) const
        {
            for( size_t i=0 ; i<x1->size() ; ++i )
                for( size_t j=0 ; j<(*x1)[0].size() ; ++j )
                    (*x1)[i][j] = m_alpha1*(*x2)[i][j] + m_alpha2*(*x3)[i][j];
        }

    };
};

#endif
