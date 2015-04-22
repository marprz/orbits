#include "AdamsIntegrator.cpp"

Vector p_Acceleration( Vector position )
{
    Vector acc = { 0, 0, 0 };
    double x = position.at(0);
    double y = position.at(1);
    double z = position.at(2);
    double r_mag = sqrt( x*x + y*y + z*z );
    double nuR3 = -kmuE/(pow(r_mag,3));
    double pp = 3/2 *kmuE*kJ2*kRE2/(pow(r_mag,5));
    double z2r2 = 5*pow(position.at(2)/r_mag,2);
    Vector c = { 1, 1, 3 };
    for( int i=0; i<3; ++i )
    {
        acc.at(i) = nuR3*position.at(i)-pp*position.at(i)*(c.at(i)-z2r2);
    }
    return acc;
}

Vector taylorTwoBody( double h, std::vector< Vector > f )
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

    Vector d3 = kmuE*( 1/r_mag3*v - 3*(r*v)*r/r_mag5);
    Vector d41 = -kmuE*( -3*r/r_mag6 * r_pr + a/r_mag3 );
    Vector d42 = 3*kmuE*( r/r_mag5*( -3*r_pr + r*a ) );
    Vector d4 = d41+d42;
    Vector ret = r + h*v + h*h*a/2 + d3*pow(h,3)/6 + d4*pow(h,4)/24;//+ d5*pow(h,5)/120;
    return ret;
}

std::vector< Vector > p_taylorTwoBodySeries( int nb, double h, std::vector< Vector > f )
{
    std::vector< Vector > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        Vector serie = taylorTwoBody( i*h, f );
        //Vector serie = taylorSerieTwoBody( i*h, f );
        ret.push_back( serie );
    }
    return ret;
}

Vector taylorTwoBodyV( double h, std::vector< Vector > f )
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

    Vector d3 = kmuE*( 1/r_mag3*v - 3*(r*v)*r/r_mag5);
    Vector d41 = -kmuE*( -3*r/r_mag6 * r_pr + a/r_mag3 );
    Vector d42 = 3*kmuE*( r/r_mag5*( -3*r_pr + r*a ) );
    Vector d4 = d41+d42;
    Vector ret = v + h*a + d3*pow(h,2)/2 + d4*pow(h,3)/6;//+ d5*pow(h,5)/120;
    return ret;
}

std::vector< Vector > p_taylorTwoBodySeriesV( int nb, double h, std::vector< Vector > f )
{
    std::vector< Vector > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        Vector serie = taylorTwoBodyV( i*h, f );
        ret.push_back( serie );
    }
    return ret;
}
void AdamsTest()
{
    AdamsIntegrator adams_integrator( p_Acceleration, p_taylorTwoBodySeries, p_taylorTwoBodySeriesV );
    double h = 15;
    int points = 1000;
    adams_integrator.Algorithm( points, h );
}

Vector p_AccSin( Vector pos )
{
    Vector ret;
    Vector::iterator it;
    for( it = pos.begin(); it != pos.end(); ++it )
        ret.push_back( -1*(*it) );
    return ret;
}

Vector taylorSin( double h, std::vector< Vector > f )
{
    Vector ret = { h-pow(h,3)/6 + pow(h,5)/120 - pow(h,7)/5040, 0, 0 };
    return ret;
}

std::vector< Vector > p_TaylorSin( int nb, double h, std::vector< Vector > f )
{
    std::vector< Vector > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        std::cout << "i*h=" << i*h << std::endl;
        double degToRad = 0.0174532925;
        Vector serie = { sin(h*i), 0, 0 };
        //Vector serie = taylorSin( i*h, f );
        ret.push_back( serie );
        Vector serie2 = { i*h-pow(i*h,3)/6+pow(i*h,5)/120-pow(i*h,7)/5040, 0, 0 };
        print( serie, "serie1" );
        print( serie2, "serie2" );
        Vector diff = serie - serie2;
        print( diff, "difference" );
    }
    return ret;
}

std::vector< Vector > p_TaylorCos( int nb, double h, std::vector< Vector > f )
{
    std::vector< Vector > ret;
    int firstIndex = -floor((nb-1)/2);
    for( int i = firstIndex; i<nb+firstIndex; ++i )
    {
        Vector serie = { cos(h*i), 0, 0 };
        ret.push_back( serie );
    }
    return ret;
}
void AdamsTestSin()
{
    Vector r0 = { 0, 0, 0 };
    Vector v0 = { 1, 0, 0 }; 
    std::cout << "AdamsTestSin: " << std::endl;
    AdamsIntegrator adams_integrator( r0, v0, p_AccSin, p_TaylorSin, p_TaylorCos );
    double h=0.05;
    int points = 500;
    std::cout << "AdamsTestSin::Algorithm: " << std::endl;
    adams_integrator.Algorithm( points, h );
}

int main()
{
    AdamsTest();
//    AdamsTestSin();
    std::cout << "TEST ADAMS INTEGRATOR" << std::endl;
    return 0;
}

