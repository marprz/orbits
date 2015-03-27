#include <iostream>

#include "Def.h"

Vector3 Acceleration( const Vector3& position )
{
    Vector3 acc = { 0, 0, 0 };
    double x = position.at(0);
    double y = position.at(1);
    double z = position.at(2);
    double r_mag = sqrt( x*x + y*y + z*z );
    double nuR3 = -kmuE/(pow(r_mag,3));
    double pp = 3/2 *kmuE*kJ2*kRE2/(pow(r_mag,5));
    double z2r2 = 5*pow(position.at(2)/r_mag,2);
    Vector3 c = { 1, 1, 3 };
    for( int i=0; i<3; ++i )
    {
        acc.at(i) = nuR3*position.at(i)-pp*position.at(i)*(c.at(i)-z2r2);
    }
    return acc;
}

int main()
{
    Vector3 p1 = { -4918.41, -14599, 21607.8 };
    Vector3 a1 = Acceleration( p1 );
    std::cout << "magnitude of acceleration: " << norm( a1 ) << std::endl;
    return 0;
}
