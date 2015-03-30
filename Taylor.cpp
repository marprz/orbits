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

template< typename T >
T taylorSerieTwoBody( double h, std::vector< T > f )
{
    T r = f.at(0);
    T v = f.at(1);
    double r_mag = norm(r);

    double product = r.at(0)*v.at(0) + r.at(1)*v.at(1) + r.at(2)*v.at(2);
    T ret = r + h*v - h*h*kmuE/(2*pow(r_mag,3)) * r + pow(h,3)*( kmuE*product*h/(2*pow(r_mag,5))*r - kmuE/(6*pow(r_mag,3))*v );
    return ret;
}

template< typename T >
std::vector< T > taylorSeriesTwoBody( int nb, double h, std::vector< T > f )
{
    std::vector< T > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        T serie = taylorSerieTwoBody( i*h, f );
        ret.push_back( serie );
    }
    return ret;
}


