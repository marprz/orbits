#include <iostream>
#include <vector>
#include <cmath>

#include "Def.h"

int fac( int a )
{
    if( a == 0 )
        return 1;
    int ret = a;
    while( --a > 1 )
        ret *= a;
    return ret;
}

// function returns values of f(x-h)
// f - vector of values: f(x), f'(x), f''(x)...
template< typename T >
T taylorSerie( double h, std::vector< T > f )
{
    T ret = f.at(0);
    for( int i=1; i<f.size(); ++i )
    {
        ret = ret + 1/fac(i)*pow(h,i)*f.at(i);
    }
    return ret;
}

template< typename T >
std::vector< T > taylorSeries( int nb, double h, std::vector< T > f )
{
    std::vector< T > ret;
    int firstIndex = -floor(nb/2);
    for( int i = firstIndex; i<(nb/2)+1; ++i )
    {
        T serie = taylorSerie( i*h, f );
        ret.push_back( serie );
    }
    return ret;
}

Vector taylorSerieTwoBody( double h, std::vector< Vector > f )
{
    Vector r = f.at(0);
    Vector v = f.at(1);
    Vector a = f.at(2);
    double r_mag = norm(r);
    double r_mag2 = r_mag*r_mag;
    double r_mag3 = r_mag2*r_mag;
    double r_mag4 = r_mag3*r_mag;
    double r_mag5 = r_mag4*r_mag;
    double r_mag6 = r_mag5*r_mag;

    double r_pr = r.at(0)*v.at(0) + r.at(1)*v.at(1) + r.at(2)*v.at(2);
    Vector vec1 = { 1, 1, 1 };
    Vector vec2 = { 2, 2, 2 };

    Vector d3 = kmuE*( 1/r_mag3*v - 3*(r*v)*r/r_mag5);
    Vector d41 = -kmuE*( -3*r/r_mag6 * r_pr + a/r_mag3 );
    Vector d42 = 3*kmuE*( r/r_mag5*( -3*r_pr + r*a ) );
    Vector d4 = d41+d42;
//    T d51 = 3/r_mag5 * v * ( 4*a - 3*r*a ) + d3/r_mag3;
//    T d52 = ( v/r_mag5 * ( vec1-5*r/r_mag2 ) )* ( r*a + v*v*(vec2-5*r/r_mag2) ) + r/r_mag5*( 5*v*a*( vec1-2*r/r_mag2 ) + r*d3 + 5*v*v*v*(2*r_mag2 - 1/r_mag2) );
//    T d5 = -kmuE*d51 + 3*kmuE*d52;
    Vector ret = r + h*v + h*h*a/2 + d3*pow(h,3)/6 + d4*pow(h,4)/24;//+ d5*pow(h,5)/120;
    return ret;
}

std::vector< Vector > taylorSeriesTwoBody( int nb, double h, std::vector< Vector > f )
{
    std::vector< Vector > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        Vector serie = taylorSerieTwoBody( i*h, f );
        ret.push_back( serie );
    }
    return ret;
}

std::vector< Vector > taylorSeriesSin( int nb, double h )
{
    std::vector< Vector > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        double h2 = h*i;
        Vector serie = { h2-pow(h2,3)/6+pow(h2,5)/120, 0, 0 };
        std::cout << "h2=" << h2 << ", taylor=" << serie.at(0) << std::endl;
        ret.push_back( serie );
    }
    return ret;
}


