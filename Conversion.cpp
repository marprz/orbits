#include <iostream>
#include <vector> 
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <tuple>
#include <queue>

#include "Def.h"

MatrixD precession( Time t, bool debug = false )
{
    MatrixD P;
    double T, z, theta, zeta;
    double JD2000_0 = 2451545; // days from Julian Date of 2000 January 1st 12h
    double OneCentury = 36525; // days of one Julian century

    // T - Julian centuries
    // funkcja JulianDate jest przetestowana :D
    double ddays = JulianDate( t.Y, t.M, t.D, t.H ) - JD2000_0;
    T = ddays/OneCentury;

    // parametry w radianach 
    z     = arcsecToRad(2306.2181*T + 1.09468*pow(T,2) + 0.018203*pow(T,3) ); 
    theta = arcsecToRad(2004.3109*T - 0.42665*pow(T,2) - 0.041833*pow(T,3) ); 
    zeta  = arcsecToRad(2306.2181*T + 0.30188*pow(T,2) + 0.017998*pow(T,3) );

    P.resize(3);
    for( int i=0; i<3; ++i )
        P.at(i).resize(3);

    P.at(0).at(0) = cos(z)*cos(theta)*cos(zeta)-sin(z)*sin(zeta);
    P.at(1).at(0) = sin(z)*cos(theta)*cos(zeta)+cos(z)*sin(zeta);
    P.at(2).at(0) = sin(theta)*cos(zeta);

    P.at(0).at(1) = -cos(z)*cos(theta)*sin(zeta)-sin(z)*cos(zeta);
    P.at(1).at(1) = -sin(z)*cos(theta)*sin(zeta)+cos(z)*cos(zeta);
    P.at(2).at(1) = -sin(theta)*sin(zeta);

    P.at(0).at(2) = -cos(z)*sin(theta);
    P.at(1).at(2) = -sin(z)*sin(theta);
    P.at(2).at(2) = cos(theta);

    if( debug )
    {
        std::cout << "Precession: T=" << T << std::endl;
        print( P, "Precession matrix" );
    }

    return P;
};

double MeanLongitudeOfAscendingLunarNode( double T );

// mean longitude of the ascending lunar node [deg]
double omegaFun( Time t )
{
    double T = JulianSince2000Century( t );
    double om = MeanLongitudeOfAscendingLunarNode( T ); // degrees
   /* double om = 125.04455501 - ( 6962890.5431*T + 7.4722*pow(T,2) + 0.007702*pow(T,3) - 0.00005939*pow(T,4) )/3600;
    double r = 1296000;
    om = 450160.280 - (5*r + 482890.539)*T + 7.455*pow(T,2) + 0.008*pow(T,3);
    if( om < 0 )
        om += r;
    om = arcsecToRad( om );*/
    om = degreeToRad( om );
    return om;
}

// obliquity of the ecliptic [degrees]
double Epsilon( Time t )
{
    double T = JulianSince2000Century(t);
    double eps = ((0.001813*T - 0.00059)*T - 46.8150)*T + 84381.448;
    eps = arcsecToDeg( eps );
    return eps;
}

  // l - Moon's mean anomaly: [deg]
double MoonMeanAnomaly( Time t )
{
    double r = 1296000;
    // obie wersje daja taki sam wynik
    double T = JulianSince2000Century(t);
    double m1 = 485866.733 + (1325*r + 715922.633)*T + 31.310*pow(T,2) + 0.064*pow(T,3);
//    double m2 = 3600*134.96340251 + 1717915923.2178*T + 31.8792*pow(T,2) + 0.051635*pow(T,3) - 0.00024470*pow(T,4);
//    std::cout << "m1 = " << m1/3600 << ", m2 = " << m2/3600 << " [deg] " << std::endl;

    // wynik jest rowny tez m1/3600 i m2/3600
    double d = JulianDate( t ) - 2451543.5;
    return 115.3654 + d*13.0649929509;
}
 
// l' - Sun's mean anomaly:
double SunMeanAnomaly( double T )
{
    double r = 1296000;
    double anomaly_arcsec = 1287099.804 + (99*r + 1292581.224)*T - 0.577*pow(T,2) - 0.012*pow(T,3);
    //return arcsecToDeg( anomaly_arcsec );
    double anomaly_deg = arcsecToDeg( anomaly_arcsec );
    return anomaly_deg;
}

// F - Moon's mean argument of latitude
double MoonMeanArgOfLatitude( double T )
{
    double r = 1296000;
    double moon_mean_latitude = 335778.877 + (1342*r + 295263.137)*T - 13.257*pow(T,2) + 0.011*pow(T,3);
    return arcsecToDeg( moon_mean_latitude );
}

// Omega - Mean longitude of the ascending lunar node:
double MeanLongitudeOfAscendingLunarNode( double T )
{
    double r = 1296000;
    double longitude = 450160.280 - (5*r + 482890.539)*T + 7.455*pow(T,2) + 0.008*pow(T,3);
    return arcsecToDeg( longitude );
}

// D - Moon's mean elongation from the Sun:
double MoonMeanElongation( double T )
{
    double r = 1296000;
    double elongation = 1072261.307 + (1236*r + 1105601.328)*T - 6.891*pow(T,2) + 0.019*pow(T,3);
    return arcsecToDeg( elongation );
}

double Mod360( double s )
{
    double ret = fmod( s, 360 );
    if( ret < 0 )
        ret += 360;
    return ret;
}

// dPsi, dEpsilon, array< double, 5 > alpha
std::tuple< double, double, std::array< double, 5 >  > NutationParameters( Time t )
{
    double T = JulianSince2000Century(t);

    double dPsi = 0, dEps = 0;

    MatrixI k;
    MatrixD A, B;
    nutationParameters( k, A, B );

    std::array< double, 5 > alpha; // degrees
    alpha.at(0) = MoonMeanAnomaly( t );
    alpha.at(1) = SunMeanAnomaly( T );
    alpha.at(2) = MoonMeanArgOfLatitude( T );
    alpha.at(3) = MoonMeanElongation( T );
    alpha.at(4) = MeanLongitudeOfAscendingLunarNode( T );

    for( int i=0; i<5; ++i )
        alpha.at(i) = degreeToRad( alpha.at(i) ); // alpha [rad]

    double eps = Epsilon( t );
    double n = k.size();
    for( int j=0; j<n; ++j )
    {
        double e1, e2=0;
        double p1, p2=0;
        e1 = B.at(j).at(0)+B.at(j).at(1)*T; // arcsec*10^(-4)
        p1 = A.at(j).at(0)+A.at(j).at(1)*T; // arcsec*10^(-4)
        for( int i=0; i<5; ++i )
        {
            p2 =+ k.at(j).at(i)*alpha.at(i); 
        }
        dEps += e1*cos(p2);
        dPsi += p1*sin(p2);

    }
    dEps *= 0.0001;
    dPsi *= 0.0001;
    dEps = arcsecToDeg( dEps );
    dPsi = arcsecToDeg( dPsi );
    return std::forward_as_tuple( dPsi, dEps, alpha );
}

// GNSS data processing - p. 177
// dpsi - temporary - TODO!
MatrixD nutation( Time t, bool debug = false )
{
    double T = JulianSince2000Century(t);

    auto parameters = NutationParameters( t );
    double dPsi = degreeToRad( std::get< 0 >( parameters ) );
    double dEpsilon = degreeToRad( std::get< 1 >( parameters ) );
    std::array< double, 5 > alpha = std::get< 2 >( parameters ); // degrees
    double epsilon = degreeToRad( Epsilon( t ) );
    double epsilonPr = epsilon + dEpsilon;

    MatrixD M1 = multiply(R3(-dPsi), R1(epsilon));
    MatrixD N = multiply(R1(-epsilonPr), M1);

    if( debug )
        print( N, "Nutation matrix" );
    return N;
}


// Earth rotation matrix
// based on article Seppanen at al: Autonomous Prediction of GPS and GLONASS
// Satellite Orbits
MatrixD rotation( Time t, bool debug = false )
{
    MatrixD G;

    // artykuł Mari
    double t_GPS = secondsSinceGPS( t ); // number of seconds elapsed since the beginning of GPS time, 6.1.1980
    double leap_seconds = leapSeconds( t );
    double JD_GPS = 2444244.5+t_GPS/86400; // [s] 
    double dUT1 = 0; // because this difference is small we neglect it
    double JD_UT1 = JD_GPS + ( dUT1 - leap_seconds )/86400;
    double JD0hUT1 = floor( JD_UT1 + 0.5 ) - 0.5; // the JD at the beginning of the current day 
    double ut1 = (JD_UT1 - JD0hUT1)*86400; // the rest of JD_UT1

    double GAST; // Greenwich Mean Sidereal Time
    double GMST; // Greenwich Apparent Sidereal Time

    double T0, T;
    T = ( JD_UT1 - 2451545 )/36525;
    T0 = ( JD0hUT1 - 2451545 )/36525;

    // GMST [seconds ]
    GMST = 24110.54841 + 1.002737909350795*ut1 + 8640184.812866*T0 + 0.093104*pow(T,2) + 6.2*pow(10,-6)*pow(T,3);
    GMST = fmod( GMST, 86400 );

    // GMST [degrees]
    GMST = 360*GMST/86400;
    // GMST [rad] 
    GMST = degreeToRad( GMST );


    double eps = Epsilon( t );
    auto parameters = NutationParameters( t );
    double dpsi = std::get< 0 >( parameters );
    double omega = omegaFun( t ); 

    dpsi = degreeToRad(-0.00222);
    double omegaE = 7.2921151*pow(10,-5); // [rad/s]

    if( debug )
        std::cout << "eps = " << eps << ", dpsi = " << dpsi << ", omega = " << omega << std::endl;

    eps = arcsecToRad( eps );
    GAST = GMST + dpsi*cos(eps);
    
    GAST = -GAST;

    G = R3( GAST );
    if( debug )
        print( G, "Rotation matrix" );
    return G;
}

struct PolarData
{
    PolarData( double at, double ax, double ay ) : t( at ), x( ax ), y( ay ) {}
    double t, x, y;
};

// polar motion matrix
MatrixD polar( Time t, bool debug = false )
{
    double x = 0; 
    double y = 0;

    std::fstream inputPolar;
    std::string fileName( "polar" );
    inputPolar.open( fileName.c_str(), std::ios::in | std::ios::out );
    int dayInYear = daysInYear( t );
    double partOfYear = static_cast< double >( dayInYear )/( t.Y%4 ? 366 : 365 );
    std::queue< PolarData > polar_xy;
    if( inputPolar.is_open() )
    {
        double xun = 0, yun =0, date = 0;
        inputPolar >> date >> x >> xun >> y >> yun;
        polar_xy.push( PolarData( date, x, y ) );
        while( inputPolar >> date >> x >> xun >> y >> yun )
        {
            if( date > (t.Y+partOfYear) )
                break;
            polar_xy.push( PolarData( date, x, y ) );
            polar_xy.back();
        }
        PolarData selected = polar_xy.back();
        inputPolar.close();
    }
    else 
    {
        std::cout << "polar IS NOT open!" << std::endl;
    }

    x = -x; y = -y;
    if( debug )
    {   
        std::cout << "-x = " << x << ", sin(x) = " << sin(degreeToRad(x)) << ", cos(x) = " << cos(degreeToRad(x)) << std::endl;
        std::cout << "-y = " << y << ", sin(y) = " << sin(degreeToRad(y)) << ", cos(y) = " << cos(degreeToRad(y)) << std::endl;
    }
    x = arcsecToRad( x );
    y = arcsecToRad( y );

    MatrixD W = multiply( R2(x), R1(y) );
    return W;
}

MatrixD conversionMatrix( Time t )
{
    MatrixD P = precession( t );
    MatrixD N = nutation( t );
    MatrixD G = rotation( t );
    MatrixD W = polar( t );

    MatrixD CEP = multiply( N, P );
    MatrixD rotated = multiply( G, CEP );
    MatrixD transformMatrix = multiply( W, rotated );

    return transformMatrix;
}

Vector3 convert( Time t, double x, double y, double z )
{
    std::vector< double > r_crs;
    r_crs.push_back( x );
    r_crs.push_back( y );
    r_crs.push_back( z );

    MatrixD M = conversionMatrix( t );
    Vector3 r_trs = multiply( M, r_crs );
    return r_trs;
}

Vector3 convert( Time t, Vector3 v )
{
    return convert( t, v.at(0), v.at(1), v.at(2) );
}

/*
int main()
{
    Time t1( 2013, 12, 10, 0, 0, 0 );
    Time t2( 2013, 12, 10, 0, 5, 0 );

    std::vector< Vector3 > pos;
    std::string filename("sat1.txt");
    std::string fileOut("out2.m");
    std::fstream input;
    std::fstream out;
    input.open( filename.c_str(), std::ios::in );
    out.open( fileOut.c_str(), std::ios::out );
    int Y = 2013;
    int M = 12;
    int D = 10; 
    int H = 0; 
    int Mi = 0; 
    int S = 0;

    out << "function A = out2\n";
    out << "A = [ ";
    //for( int i=0; i<2; ++i )
    for( int i=0; i<24*12; ++i )
    {
        double x, y, z; 
        std::string str;
        input >> x >> y >> z >> str;
        std::cout << "Coordinates: " << x << ", " << y << ", " << z << std::endl;
        Time t( Y, M, D, H, Mi, S );
        Vector3 v = convert( t, x, y, z, fmod(H,12)+Mi/60 );
        pos.push_back( v );
        if( Mi < 55 )
        {
            Mi += 5;
        }
        else
        {
            Mi = 0;
            ++H;
        }
        out << (int)v.at(0) << " " << (int)v.at(1) << " " << (int)v.at(2) << ";\n";
    }
    out << "];";
    out.close();
    input.close();


/*    std::cout << "First conversion: " << std::endl;
    convert( t1, -15210.576086, 2042.480101, 21649.168938 );
    std::cout << "Second conversion: " << std::endl;
    convert( t2, 14245.419967, 20412.183231, -9713.789057 );*/
/*    return 0;
}*/
