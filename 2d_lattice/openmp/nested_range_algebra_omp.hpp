/* nested range algebra */

#ifndef NESTED_OMP_ALGEBRA_HPP
#define NESTED_OMP_ALGEBRA_HPP

template< class InnerAlgebra >
struct nested_omp_algebra
{
    
    nested_omp_algebra( const int block_size = 1 ) 
        : m_inner_algebra() , m_block_size( block_size )
    { }

    template< class S1 , class S2 , class S3 , class Op >
    void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
    {
#pragma omp for schedule( dynamic , m_block_size )
        for( size_t i=0 ; i<boost::size(s1) ; ++i )
            m_inner_algebra.for_each3( s1[i] , s2[i] , s3[i] , op );
    }


private:
    InnerAlgebra m_inner_algebra;
    const int m_block_size;
};

#endif
