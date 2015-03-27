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
double taylorSerie( double h, std::vector< double > f )
{
    double ret = f.at(0);
    for( int i=1; i<f.size(); ++i )
    {
        ret += f.at(i)*pow(h,i)/fac(i);
    }
    return ret;
}

std::vector< double > taylorSeries( int nb, double h, std::vector< double > f )
{
    std::vector< double > ret;
    int firstIndex = 0-floor(nb/2);
    for( int i = firstIndex; i<(nb/2)+1; ++i )
    {
        double serie = taylorSerie( i*h, f );
        ret.push_back( serie );
    }
    return ret;
}


