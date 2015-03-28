#include "Integration.h"
#include "Acceleration.h"
#include "Conversion.cpp"

void SimulationFromZero()
{
    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    Vector3 r0 = { rp, 0, 0 };
    Vector3 v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };

    // vectors of positions, velocities
    std::vector< Vector3 > rs;
    std::vector< Vector3 > vs;

    rs.push_back( r0 );
    vs.push_back( v0 );

    double h = 60; // step size in seconds
    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector3 > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector3 > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positions");

}

// r - position
// v - velocity 
// a - acceleration
void SimulationNGA()
{
    // vectors of positions, velocities, accelerations
    std::vector< Vector3 > rs;
    std::vector< Vector3 > vs;
    std::vector< Vector3 > as;

    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h
//    double Vp = sqrt( muE*(2/rp-1/sa)); // speed at perigee

    Time t( 2013, 12, 10, 0, 0, 0 );
    Vector3 r0_TRS = { -15210.576086, 2042.480101, 21649.168938 };
    Vector3 r0 = convert( t, r0_TRS );
    Vector3 v0_TRS = { -7773.910720, -26449.272652, -2851.405296 };
    Vector3 v0 = 0.001*v0_TRS ;
    //Vector3 v0 = convert( t, v0_TRS );

    double r = norm(r0);
    std::cout << "Vp = " << sqrt( kmuE*(2/r-1/sa) ) << std::endl;
    std::cout << "v0 = " << v0.at(0) << " " << v0.at(1) << " " << v0.at(2) << std::endl;
    rs.push_back( r0 );
    vs.push_back( v0 );
    as.push_back( Acceleration( r0 ) );

    double h = 60; // step size: 5 minutes
    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector3 > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector3 > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positions");
}

void SimulationSin()
{
    std::vector< Vector3 > rs;
    std::vector< Vector3 > vs;
  
    Vector3 r0 = { 0 };
    Vector3 v0 = { 1 };
    rs.push_back( r0 );
    vs.push_back( v0 );
    double h = 0.02;
    for( int i=0; i<300; ++i )
    {
        std::vector< Vector3 > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }
    std::cout << "wektor pozycji dla sinusa (h=0.1): ";
    for( int i=0; i<rs.size(); ++i )
    {
        std::cout << rs.at(i).at(0) << " (" << sin(i*h) << "), ";
    }
    std::cout << std::endl;
}

int main()
{
    SimulationFromZero();
//    SimulationNGA(); 
//    SimulationSin();
    return 0;
}
