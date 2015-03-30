#include <list>
#include "Integration.h"
#include "Acceleration.h"
#include "Conversion.cpp"
#include "Taylor.cpp"

typedef std::vector< std::vector< Vector > > InitialVectors;

InitialVectors InitializeVectorsGJ( int nb, double h )
{
    InitialVectors init;
    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };
    Vector a0 = Acceleration( r0 );

    std::vector< Vector > fr = { r0, v0, a0 };
    std::vector< Vector > fv = { v0, a0 };
    // TWO-BODY TAYLOR:
    std::vector< Vector > rs = taylorSeriesTwoBody( 10, h, fr ); // potrzeba 9 pozycji do algorytmu. 10-ta jest do obliczenia predkosci jako pochodnych pozycji
    std::vector< Vector > vs;
    std::vector< Vector > as;
    std::vector< Vector >::iterator it;
    for( int i=0; i<rs.size()-1; ++i )
    {
        vs.push_back( 1/h*(rs.at(i+1)-rs.at(i)) );
    }
    rs.pop_back(); // ostatni element byl potrzebny tylko do liczenia predkosci
    for( it = rs.begin(); it != rs.end(); ++it )
    {
        as.push_back( Acceleration( *it ) );
    }
    init.push_back( rs );
    init.push_back( vs );
    init.push_back( as );
    return init;
}

// equation 75 s0 = C1'
Vector Calcs0( const Vector& v0, const MatrixD& b, const std::vector< Vector >& as, const double& h )
{
    Vector ret;
    Vector sum = { 0, 0, 0 };
    for( int i=0; i<as.size(); ++i )
        sum = sum + b.at(4).at(i)*as.at(i);
    ret = v0/h-sum;
    print( ret, "s0" );
    return ret;
}

// equation 86
Vector CalcS0( const Vector& v0, const MatrixD& a, const std::vector< Vector >&as, const double& h )
{
    Vector ret;
    Vector sum = { 0, 0, 0 };
    for( int i=0; i<as.size(); ++i )
        sum = sum + a.at(4).at(i)*as.at(i);
    ret = v0/pow(h,2)-sum;
    print( ret, "s0" );
    return ret; 
}
// equation 75
void Calcsn( std::list< Vector >& s, const std::vector< Vector >& a, const int& index )
{
    if( index == 0 )
        std::cout << "ERROR! Calcsn - index should not be 0! " << std::endl;
    if( index < 0 )
    {
        // example: index=-1 => a(n) = 3, a(n+1) = 4; 3=index+4, 4=index+5
        Vector sn = s.front() - ( a.at(index+5) + a.at(index+4) )/2;
        s.push_front( sn );
    } 
    else
    {
        Vector sn = s.back() + ( a.at(index+4-1) + a.at(index+5-1) )/2;
        s.push_back( sn );
    }
}

// equation 86
void CalcSn( const std::vector< Vector >& s, std::list< Vector >& S, const std::vector< Vector >& a, const int& index )
{
    if( index == 0 )
        std::cout << "ERROR! CalcSn - index should not be 0! " << std::endl;
    if( index < 0 )
    {
        // example: index=-1 => a(n) = 3, a(n+1) = 4; 3=index+4, 4=index+5
        Vector Sn = S.front() - s.at(index+5) +(a.at(index+5))/2;
        S.push_front( Sn );
    } 
    else
    {
        Vector Sn = S.back() + s.at(index+4) + a.at(index+4)/2;
        S.push_back( Sn );
    }
}

// table 5
std::vector< Vector > CalcSumb( const MatrixD& b, std::vector< Vector > as )
{
    std::vector< Vector > ret;
    for( int i=0; i<9; ++i ) // n=-4,...,4
    {
        Vector sum = { 0, 0, 0 };
        for( int k=0; k<9; ++k )
        {
            sum = sum + b.at(i).at(k)*as.at(k);
        }
        ret.push_back( sum );
    }
    return ret;
}

// table 6
std::vector< Vector > CalcSuma( const MatrixD& a, std::vector< Vector > as )
{
    std::vector< Vector > ret;
    for( int i=0; i<9; ++i ) // n=-4,...,4
    {
        Vector sum = { 0, 0, 0 };
        for( int k=0; k<9; ++k )
        {
            sum = sum + a.at(i).at(k)*as.at(k);
        }
        ret.push_back( sum );
    }
    return ret;
}



void SimulationFromZeroGaussJackson()
{    int n = 9;
    double h = 60; // step size in seconds
    InitialVectors init = InitializeVectorsGJ( n, h );
    std::vector< Vector > rs = init.at(0);
    std::vector< Vector > vs = init.at(1);
    std::vector< Vector > as = init.at(2);

    std::cout << "startup - r, v, a sizes: " << rs.size() << " " << vs.size() << " " << as.size() << std::endl;
    print( rs, "two-body positions");
    print( vs, "two-body velocities");
    print( as, "two-body accelerations");
    std::vector< Vector >::iterator it;
    std::cout << "positions: ";
    for( it = rs.begin(); it != rs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << "velocities: ";
    for( it = vs.begin(); it != vs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << "accelerations: ";
    for( it = as.begin(); it != as.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;

    MatrixD a, b;
    GaussJacksonCoefficients( a, b );

    bool is_acceleration_converged = false;
    while( !is_acceleration_converged )
    {
        std::list< Vector > snList, SnList;
        std::vector< Vector > sum_bn, sum_an;

        // sn:
        Vector s0 = Calcs0( vs.at(4), b, as, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
        snList.push_back( s0 );
        for( int i=1; i<5; ++i )
        {
            //Calcsn( std::list< Vector >& s, const std::vector< Vector >& a, const int& index )
            Calcsn( snList, as, -i ); 
            Calcsn( snList, as, i ); 
        }

        std::vector< Vector > sn, Sn; // TODO powinna byc lepsza struktura, zamiast tworzenia nowego wektora z listy
        for( int j=0; j<9; ++j )
        {
            sn.push_back( snList.front() );
            snList.pop_front();
        }

        // Sn:
        Vector S0 = CalcS0( vs.at(4), a, as, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
        SnList.push_back( S0 );
        for( int i=1; i<5; ++i )
        {
            CalcSn( sn, SnList, as, -i ); // CalcSn( const std::vector< Vector >& s, std::list< Vector >& S, const std::vector< Vector >& a, const int& index )
            CalcSn( sn, SnList, as, i ); 
        }
        
        for( int j=0; j<9; ++j )
        {
            Sn.push_back( SnList.front() );
            SnList.pop_front();

        }
        print( sn, "sn" );
        print( Sn, "Sn" );

        // step 3b  iii.
        sum_bn = CalcSumb( b, as ); 
        sum_an = CalcSuma( a, as );

        // step 3b iv 
        std::vector< Vector > corrected_vs; // eq. 74
        std::vector< Vector > corrected_rs; // eq. 87
        for( int i=0; i<9; ++i )
        {   
            if( i == 4 )
            {
                corrected_vs.push_back( vs.at(4) );
                corrected_rs.push_back( rs.at(4) );
            }
            else 
            {
                corrected_vs.push_back( h*(sn.at(i)+sum_bn.at(i)) );
                corrected_rs.push_back( h*h*(Sn.at(i)+sum_an.at(i)) );
            }
        }
        print( rs, "two-body positions");
        print( corrected_rs, "corrected_rs" );
        print( vs, "two-body velocities");
       // print( as, "two-body accelerations");

        print( corrected_vs, "corrected_vs" );

        is_acceleration_converged = true; // TODO temporary
    }




/*    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positions");
*/
}

void SimulationFromZeroNystrom()
{
    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };

    // vectors of positions, velocities
    std::vector< Vector > rs;
    std::vector< Vector > vs;

    rs.push_back( r0 );
    vs.push_back( v0 );

    double h = 60; // step size in seconds
    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positions");

}

// r - position
// v - velocity 
// a - acceleration
void SimulationNGA() // nie dziala
{
    // vectors of positions, velocities, accelerations
    std::vector< Vector > rs;
    std::vector< Vector > vs;
    std::vector< Vector > as;

    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h
//    double Vp = sqrt( muE*(2/rp-1/sa)); // speed at perigee

    Time t( 2013, 12, 10, 0, 0, 0 );
    Vector r0_TRS = { -15210.576086, 2042.480101, 21649.168938 };
    Vector r0 = convert( t, r0_TRS );
    Vector v0_TRS = { -7773.910720, -26449.272652, -2851.405296 };
    Vector v0 = 0.001*v0_TRS ;
    //Vector v0 = convert( t, v0_TRS );

    double r = norm(r0);
    std::cout << "Vp = " << sqrt( kmuE*(2/r-1/sa) ) << std::endl;
    std::cout << "v0 = " << v0.at(0) << " " << v0.at(1) << " " << v0.at(2) << std::endl;
    rs.push_back( r0 );
    vs.push_back( v0 );
    as.push_back( Acceleration( r0 ) );

    double h = 60; // step size: 5 minutes
    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positions");
}

void SimulationSin()
{
    std::vector< Vector > rs;
    std::vector< Vector > vs;
  
    Vector r0 = { 0 };
    Vector v0 = { 1 };
    rs.push_back( r0 );
    vs.push_back( v0 );
    double h = 0.02;
    for( int i=0; i<300; ++i )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
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
    SimulationFromZeroGaussJackson();
//    SimulationFromZeroNystrom();
//    SimulationNGA();  // nie dziala - V0 nie jest znane
//    SimulationSin();
    return 0;
}
