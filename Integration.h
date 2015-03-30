#include "Def.h"
#include "Acceleration.h"

std::vector< Vector > GaussJackson( const Vector& r0, const Vector& v0, const double&h )
{

}

// return vector: r1, v1
std::vector< Vector > OdeNystrom( const Vector& r0, const Vector& v0, const double& h )
{
    Vector ALPHA;
    std::vector< Vector > A;
    Vector AL;

    HairerParameters( ALPHA, A, AL );

    Vector tempV = ALPHA.at(0)*h*v0;
    tempV = r0+tempV;
//    Vector temp = { -1*r0 };
    Vector temp = Acceleration( tempV );

    const int s = 12; // from Hairer
    std::array< Vector, s > k;
    k.at(0) = temp;
    for( int i=1; i<s-1; ++i )
    {
        Vector sum = { 0, 0, 0 };
        for( int j=0; j<i-1; ++j )
        {
            sum = sum + (A.at(i-1).at(j))*(k.at(j));
        }
        Vector par = pow(h,2)*sum + r0+ALPHA.at(i)*h*v0;
        k.at(i) = Acceleration( par );
//        k.at(i) = -1*par;
    }

    Vector sum2 = { 0, 0, 0 };
    for( int i=0; i<s-1; ++i )
    {
        sum2 = sum2 + A.at(s-1).at(i)*k.at(i);
    }
    Vector r1 = h*h*sum2 + r0+h*v0;

    Vector sum3 = { 0, 0, 0 };
    for( int i=0; i<s-1; ++i )
    {
        sum3 = sum3 + AL.at(i)*k.at(i);
    }
    Vector v1 = h*sum3+v0;

    std::vector< Vector > ret;
    ret.push_back( r1 );
    ret.push_back( v1 );

    return ret;
    
}
