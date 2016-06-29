#include "NystromIntegrator.h"

RungeKutta::RungeKutta()
{
    ra = kRE+20394; // apogee [km]
    rp = kRE+19970; // perigee [km]
    sa = (rp+ra)/2;  // semi-major axis
    T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    std::cout << "perigee: " << rp << std::endl;
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };

    Startup( r0, v0 );
}

void NystromIntegrator::Startup( const Vector& position, const Vector& velocity )
{
    rs.push_back( position );
    vs.push_back( velocity );
    as.push_back( Acceleration( position ) );
}

Vector NystromIntegrator::Acceleration( const Vector& position )
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

void NystromIntegrator::Algorithm( int points, double h )
{
    for( int i=0; i<points; ++i )
    {
        Vector r0 = rs.back();
        Vector v0 = vs.back();
       rs.push_back( r1 );
        vs.push_back( v1 );
        as.push_back( Acceleration( r1 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positionsNystrom");

}

