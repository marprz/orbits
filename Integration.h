#include "Def.h"
#include "Acceleration.h"

std::vector< Vector3 > GaussJackson( const Vector3& r0, const Vector3& v0, const double&h )
{

}

// return vector: r1, v1
std::vector< Vector3 > OdeNystrom( const Vector3& r0, const Vector3& v0, const double& h )
{
    Vector3 ALPHA;
    std::vector< Vector3 > A;
    Vector3 AL;

    HairerParameters( ALPHA, A, AL );

    Vector3 tempV = ALPHA.at(0)*h*v0;
    tempV = r0+tempV;
//    Vector3 temp = { -1*r0 };
    Vector3 temp = Acceleration( tempV );

    const int s = 12; // from Hairer
    std::array< Vector3, s > k;
    k.at(0) = temp;
    for( int i=1; i<s-1; ++i )
    {
        Vector3 sum = { 0, 0, 0 };
        for( int j=0; j<i-1; ++j )
        {
            sum = sum + (A.at(i-1).at(j))*(k.at(j));
        }
        Vector3 par = pow(h,2)*sum + r0+ALPHA.at(i)*h*v0;
        k.at(i) = Acceleration( par );
//        k.at(i) = -1*par;
    }

    Vector3 sum2 = { 0, 0, 0 };
    for( int i=0; i<s-1; ++i )
    {
        sum2 = sum2 + A.at(s-1).at(i)*k.at(i);
    }
    Vector3 r1 = h*h*sum2 + r0+h*v0;

    Vector3 sum3 = { 0, 0, 0 };
    for( int i=0; i<s-1; ++i )
    {
        sum3 = sum3 + AL.at(i)*k.at(i);
    }
    Vector3 v1 = h*sum3+v0;

    std::vector< Vector3 > ret;
    ret.push_back( r1 );
    ret.push_back( v1 );

    return ret;
    
}
