#ifndef DEF_H
#define DEF_H

#include <array>
#include <vector>
#include <cmath>
#include <iomanip>

typedef std::array< std::array< double, 3 >, 3 > Matrix33;
typedef std::vector< std::vector< double > > MatrixD;
typedef std::vector< std::vector< int > > MatrixI;
typedef std::vector< double > Vector3;

double kmuE = 398600.4418;
double kRE = 6738;
double kRE2 = kRE*kRE;
double kJ2 = 0.00108263;

MatrixD createMatrix( std::vector< double > v )
{
    MatrixD M; 
    Vector3 v1;
    v1.push_back( v.at(0) );
    v1.push_back( v.at(1) );
    v1.push_back( v.at(2) );
    Vector3 v2;
    v2.push_back( v.at(3) );
    v2.push_back( v.at(4) );
    v2.push_back( v.at(5) );
    Vector3 v3;
    v3.push_back( v.at(6) );
    v3.push_back( v.at(7) );
    v3.push_back( v.at(8) );
    M.push_back( v1 );
    M.push_back( v2 );
    M.push_back( v3 );

    return M;
}

void print( MatrixD, std::string );

MatrixD InverseMatrix( const MatrixD& M )
{
    double a = M.at(0).at(0);
    double b = M.at(0).at(1);
    double c = M.at(0).at(2);
    double d = M.at(1).at(0);
    double e = M.at(1).at(1);
    double f = M.at(1).at(2);
    double g = M.at(2).at(0);
    double h = M.at(2).at(1);
    double i = M.at(2).at(2);

    double A = e*i-f*h;
    double B = -(d*i-f*g);
    double C = d*h-e*g;
    double D = -(b*i-c*h);
    double E = a*i-c*g;
    double F = -(a*h-b*g);
    double G = b*f-c*e;
    double H = -(a*f-c*d);
    double I = a*e-b*d;

    double detM = a*A+b*B+c*C;
    detM = a*e*i + b*f*g + d*d*h - c*d*g-f*h*a-i*b*d;
    if( detM == 0 )
        std::cout << "ERROR: INVERSE MATRIX DOES NOT EXIST" << std::endl;

    std::cout << "det = " << detM << std::endl;
    Vector3 v1 = { A/detM, D/detM, G/detM };
    Vector3 v2 = { B/detM, E/detM, H/detM };
    Vector3 v3 = { C/detM, F/detM, I/detM };

    MatrixD inverse = { v1, v2, v3 };
    return inverse;
}

void Transpose( MatrixD& M )
{
    for( int i=0; i<M.size(); ++i )
    {
        for( int j=i; j<M.size(); ++j )
        {
            double temp = M.at(i).at(j); 
            M.at(i).at(j) = M.at(j).at(i);
            M.at(j).at(i) = temp;
        }
    }
}

void print( Vector3 V, std::string name = "Vector" )
{
    std::cout << std::endl << name << ": ";
//    double r = sqrt( pow(V.at(0),2) + pow(V.at(1),2) + pow(V.at(2),2) );
    for( int i=0; i<V.size(); ++i )
        std::cout << V.at(i) << " ";
    std::cout << std::endl;
}

double degreeToRad( double deg )
{
    double degToRad = 0.0174532925;
    return deg*degToRad;
}

double radToDegree( double rad )
{
    double degToRad = 0.0174532925;
    return rad/degToRad;
}

// arcsecond symbol: ''
// http://en.wikipedia.org/wiki/Minute_of_arc
double arcsecToDeg( double s )
{
    return s/3600;
}

double arcsecToRad( double s )
{
    return degreeToRad( arcsecToDeg( s ) );
}

double arcsecToDegMod30( double s )
{
    double deg = arcsecToDeg( s );
    deg = fmod( deg, 360 );
    if( deg < 360 )
        deg += 360;
    return deg;
}

// TODO przerobic drukowanie na template
void print( MatrixD M, std::string name = "Matrix" )
{
    std::cout << std::endl << name << ": " << std::endl;
    std::cout << M.at(0).at(0) << " " << M.at(0).at(1) << " " << M.at(0).at(2) << std::endl;
    std::cout << M.at(1).at(0) << " " << M.at(1).at(1) << " " << M.at(1).at(2) << std::endl;
    std::cout << M.at(2).at(0) << " " << M.at(2).at(1) << " " << M.at(2).at(2) << std::endl;
}

void print( Matrix33 M, std::string name = "Matrix" )
{
    std::cout << std::endl << name << ": " << std::endl;
    std::cout << M.at(0).at(0) << " " << M.at(0).at(1) << " " << M.at(0).at(2) << std::endl;
    std::cout << M.at(1).at(0) << " " << M.at(1).at(1) << " " << M.at(1).at(2) << std::endl;
    std::cout << M.at(2).at(0) << " " << M.at(2).at(1) << " " << M.at(2).at(2) << std::endl;
}

// @param m - angle in radians
MatrixD R1( double m )
{
    MatrixD R;

    Vector3 v1;
    v1.push_back( 1 );
    v1.push_back( 0 );
    v1.push_back( 0 );
   
    Vector3 v2;
    v2.push_back( 0 );
    v2.push_back( cos(m) );
    v2.push_back( sin(m) );

    Vector3 v3;
    v3.push_back( 0 );
    v3.push_back( -sin(m) );
    v3.push_back( cos(m) );

    //      | v1.at(0) v1.at(1) v1.at(2) |
    // R =  | v2.at(0) v2.at(1) v2.at(2) |
    //      | v3.at(0) v3.at(1) v3.at(2) |
    R.push_back( v1 );
    R.push_back( v2 );
    R.push_back( v3 );

    return R;
}

MatrixD R2( double m )
{
    MatrixD R;

    Vector3 v1;
    v1.push_back( cos(m) );
    v1.push_back( 0 );
    v1.push_back( -sin(m) );
   
    Vector3 v2;
    v2.push_back( 0 );
    v2.push_back( 1 );
    v2.push_back( 0 );

    Vector3 v3;
    v3.push_back( sin(m) );
    v3.push_back( 0 );
    v3.push_back( cos(m) );

    R.push_back( v1 );
    R.push_back( v2 );
    R.push_back( v3 );

    return R;
}

MatrixD R3( double m )
{
    MatrixD R;

    Vector3 v1;
    v1.push_back( cos(m) );
    v1.push_back( sin(m) );
    v1.push_back( 0 );
   
    Vector3 v2;
    v2.push_back( -sin(m) );
    v2.push_back( cos(m) );
    v2.push_back( 0 );

    Vector3 v3;
    v3.push_back( 0 );
    v3.push_back( 0 );
    v3.push_back( 1 );

    R.push_back( v1 );
    R.push_back( v2 );
    R.push_back( v3 );

    return R;
}
/*
Matrix33 R3( double m )
{
    Matrix33 R;
    R.at(0).at(0) = cos( m );
    R.at(1).at(0) = -sin( m );
    R.at(2).at(0) = 0;
    
    R.at(0).at(0) = sin( m );
    R.at(1).at(0) = cos( m );
    R.at(2).at(0) = 0;

    R.at(0).at(0) = 0;
    R.at(1).at(0) = 0;
    R.at(2).at(0) = 1;

    std::cout << "testowa macierz dla m=" << m << std::endl;
    std::cout << "cos(m) = " << cos( m ) << std::endl;
    std::cout << "sin(m) = " << sin( m ) << std::endl;
    for( int i=0; i<3; ++i )
        std::cout << R.at(i).at(0) << " " << R.at(i).at(1) << " " << R.at(i).at(2) << std::endl;
    print( R, "testowa R3" );
    return R;
}*/

MatrixD multiply( MatrixD M, MatrixD N )
{
    MatrixD ret;
    for( int i=0; i<3; ++i )
    {
        Vector3 v;
        for( int j=0; j<3; ++j )
        {
            double s = 0;
            for( int k=0; k<3; ++k )
                s += M.at(i).at(k)*N.at(k).at(j);
            v.push_back( s );
        }
        ret.push_back( v );
    }
    return ret;
}

Vector3 multiply( MatrixD M, Vector3 v )
{
    Vector3 ret;
    for( int i=0; i<3; ++i )
    {
        double s = 0;
        for( int j=0; j<3; ++j )
            s += v.at(j)*M.at(i).at(j);
        ret.push_back( s );
    }
    return ret;
}

struct Time
{
    Time( int y, int m, int d, int h, int mt, int s ) : Y(y), M(m), D(d), H(h), Mt(mt), S(s) {}
    Time() : Y(0), M(0), D(0), H(0), Mt(0), S(0) {}
    void print( std::string a_name );
    int Y, M, D, H, Mt, S; 
};

void Time::print( std::string a_name = "" )
{
    std::cout << "Date " << a_name << ": " << D << "." << M << "." << Y << ", " << H << ":" << Mt << ":" << S << std::endl;
}

// works only after 1980 and before 2100 
int daysInYear( Time t, bool leapYear = false )
{
    leapYear = !(t.Y%4);
    int ret = t.D;
    if( t.M > 1 ) 
        ret += 31;
    if( t.M > 2 )
        leapYear ? ret+=29 : ret+=28; 
    if( t.M > 3 )
        ret += 31;
    if( t.M > 4 )
        ret += 30;
    if( t.M > 5 )
        ret += 31;
    if( t.M > 6 )
        ret += 30;
    if( t.M > 7 )
        ret += 31;
    if( t.M > 8 )
        ret += 31;
    if( t.M > 9 )
        ret += 30;
    if( t.M > 10 )
        ret += 31;
    if( t.M > 11 )
        ret += 30;
    return ret;
}

int leapYearsSince1980( int y )
{
    int ret = 0;
    y -= 1980;
    ret = (y-1)/4; // TODO sprawdzic, czy na pewno trzeba odjac '1'
    return ret;
}

// number of seconds since 6 January 1980, 
int secondsSinceGPS( Time t )
{
    int days = 0;
    int seconds = 0;
    days += daysInYear( t, !(t.Y%4) );
    days -= 5;
    days += (t.Y-1980)*365;
    days += leapYearsSince1980( t.Y );
    seconds = days*24*60*60;
    seconds += t.H*60*60;
    seconds += t.Mt*60;
    seconds += t.S;
    return seconds;
}

int leapSeconds( Time t )
{
    if( t.Y <= 1972 )
    {
        if( t.M < 7 )
            return 0;
        else return 1;
    }
    if( t.Y > 1972 && t.Y < 1980 )
        return (2+t.Y-1973);
    if( t.Y == 1980 )
        return 9;
    if( t.Y < 1984 )
        if( t.M > 6 )
            return (9+t.Y-1980);
        else
            return (9+t.Y-1981);
    if( t.Y == 1984 )
        return 12;
    if( t.Y == 1985 )
        if( t.M <= 6 )
            return 12;
        else
            return 13;
    if( t.Y == 1986 || t.Y == 1987 )
        return 13;
    if( t.Y == 1988 || t.Y == 1989 )
        return 14;
    if( t.Y == 1990 )
        return 15;
    if( t.Y == 1991 || ( t.Y == 1992 && t.M <= 6 ) )
        return 16;
    if( t.Y == 1992 || ( t.Y == 1993 && t.M <= 6 ) )
        return 17;
    if( t.Y == 1993 || ( t.Y == 1994 && t.M <= 6 ) )
        return 18;
    if( t.Y == 1994 || t.Y == 1995 )
        return 19;
    if( t.Y == 1996 || ( t.Y == 1997 && t.M <= 6 ) )
        return 20;
    if( t.Y == 1997 || t.Y == 1998 )
        return 21;
    if( t.Y > 1998 && t.Y <= 2005 )
        return 22;
    if( t.Y > 2005 && t.Y <= 2008 )
        return 23;
    if( t.Y > 2008 && ( t.Y <= 2012 && t.M <= 6 ) ) 
        return 24;
    if( t.Y >= 2012 )
        return 25;

    std::cout << "Blad w leapSeconds" << std::endl;
    return 0;
        

}
// conversion calendar date to Julian date
// day, month, year, hour
double JulianDate( int Y, int M, int D, int H )
{
    double JD;
    int y, m; 
    if( M > 2 )
    {
        y = Y;
        m = M;
    }
    else
    {
        y = Y-1;
        m = M+12;
    }

    double p1 = 365.25*y;
    double p2 = (m+1)*30.6001;

    JD = static_cast< int >(365.25*y) + static_cast< int >( 30.6001*(m+1) ) + D + (double)H/24 + 1720981.5;
    //JD = p1 + p2 + D + (double)H/24+1720981.5;

    return JD;
}

double JulianDate( Time t )
{
    return JulianDate( t.Y, t.M, t.D, t.H );
}

double JulianSince2000Century( Time t )
{
    int OneCentury = 36525;
    double jd1 = JulianDate( t.Y, t.M, t.D, t.H );
    jd1 /= OneCentury;
    Time t2( 2000, 1, 1, 12, 0, 0 );
    double jd2 = JulianDate( t2.Y, t2.M, t2.D, t2.H );
    jd2 /= OneCentury;

    return jd1-jd2;
    //return jd2-jd1;
}

void inverseJulian( double jd )
{
    int a = static_cast< int >( jd+0.5 );
    int b = a+1537;
    int c = static_cast< int >( (b-122.1)/365.25 );
    int d = static_cast< int >( 365.25*c );
    int e = static_cast< int >( (b-d)/30.6001 );

    int D = b-d-static_cast< int >( 30.6001*e )+( jd+0.5 ) - static_cast< int >( jd+0.5 );
    int M = e-1-12*static_cast< int >( e/14 );
    int Y = c-4715-static_cast< int >( (7+M)/10 );
    double H = (jd+0.5-static_cast< int >(jd+0.5))*24;

    std::cout << "inverse Julian: " << Y << "." << M << "." << D << ", H:" << H << std::endl;
}

Vector3 operator- ( const Vector3& v1, const Vector3& v2 )
{
    Vector3 ret;
    for( int i=0; i<3; ++i )
        ret.push_back( v1.at(i)-v2.at(i) );
    return ret;
}

double norm( const Vector3& v )
{
    double ret = 0;
    for( int i=0; i<v.size(); ++i )
        ret += (v.at(i))*(v.at(i));
    ret = sqrt( ret );
    return ret;
}

// vector: <x,y,z>
// parameter: vector < R, Lat, Lon > 
std::vector< double > LatLonToXYZ( std::vector< double > RLatLon )
{
    std::vector< double > ret; 
    double x,y,z;
    double r, lat, lon;
    r = RLatLon.at(0);
    lat = degreeToRad( RLatLon.at(1) );
    lon = degreeToRad( RLatLon.at(2) );
    x = r*cos(lat)*cos(lon);
    y = r*cos(lat)*sin(lon);
    z = r*sin(lat);
    ret.push_back( x );
    ret.push_back( y );
    ret.push_back( z );
    return ret;
};

// GNSS data processing - p.181
void nutationParameters( MatrixI& k, MatrixD& A, MatrixD& B )
{
    k.push_back( std::vector< int >{ 0, 0, 0, 0, 1 } ); //1
    k.push_back( std::vector< int >{ 0, 0, 2, -2, 2} );
    k.push_back( std::vector< int >{ 0, 0, 2, 0, 2 } );
    k.push_back( std::vector< int >{ 0, 0, 0, 0, 2 } );
    k.push_back( std::vector< int >{ 0, -1, 0, 0, 0} ); //5
    k.push_back( std::vector< int >{ 1, 0, 0, 0, 0} );
    k.push_back( std::vector< int >{ 0, 1, 2, -2, 1} );
    k.push_back( std::vector< int >{ 0, 0, 2, 0, 1} );
    k.push_back( std::vector< int >{ 1, 0, 2, 0, 2} );
    k.push_back( std::vector< int >{ 0, -1, 2, -2, 2} ); // 10
    k.push_back( std::vector< int >{ -1, 0, 0, 2, 0} );
    k.push_back( std::vector< int >{ 0, 0, 2, -2, 1} );
    k.push_back( std::vector< int >{ -1, 0, 2, 0, 2} );
    k.push_back( std::vector< int >{ 1, 0, 0, 0, 1} );
    k.push_back( std::vector< int >{ 0, 0, 0, 2, 0} ); // 15
    k.push_back( std::vector< int >{ -1, 0, 2, 2, 2} );
    k.push_back( std::vector< int >{ -1, 0, 0, 0, 1} );
    k.push_back( std::vector< int >{ 1, 0, 2, 0, 1} );
    k.push_back( std::vector< int >{ -2, 0, 0, 2, 0} );
    k.push_back( std::vector< int >{ -2, 0, 2, 0, 1} ); // 20 
    k.push_back( std::vector< int >{ 0, 0, 2, 2, 2} );
    k.push_back( std::vector< int >{ 2, 0, 2, 0, 2} );
    k.push_back( std::vector< int >{ 2, 0, 0, 0, 0} );
    k.push_back( std::vector< int >{ 1, 0, 2, -2, 2} );
    k.push_back( std::vector< int >{ 0, 0, 2, 0, 0} ); // 25
    k.push_back( std::vector< int >{ 0, 0, 2, -2, 2} );
    k.push_back( std::vector< int >{ -1, 0, 2, 0, 1} );
    k.push_back( std::vector< int >{ 0, 2, 0, 0, 0} );
    k.push_back( std::vector< int >{ 0, 2, 2, -2, 2} );
    k.push_back( std::vector< int >{ -1, 0, 0, 2, 1} ); // 30
    k.push_back( std::vector< int >{ 0, 1, 0, 0, 1} );
    k.push_back( std::vector< int >{ 1, 0, 0, -2, 1} );
    k.push_back( std::vector< int >{ 0, -1, 0, 0, 1} );
    k.push_back( std::vector< int >{ 2, 0, -2, 0, 0} );
    k.push_back( std::vector< int >{ -1, 0, 2, 2, 1} ); // 35
    k.push_back( std::vector< int >{ 1, 0, 2, 2, 2} );
    k.push_back( std::vector< int >{ 0, -1, 2, 0, 2} );
    k.push_back( std::vector< int >{ 0, -1, 2, 0, 2} );
    k.push_back( std::vector< int >{ 1, 1, 0, -2, 0} );
    k.push_back( std::vector< int >{ 0, 1, 2, 0, 2} ); // 40
    k.push_back( std::vector< int >{ -2, 0, 0, 2, 1 } ); 
    k.push_back( std::vector< int >{ 0, 0, 0, 2, 1 } );
    k.push_back( std::vector< int >{ 2, 0, 2, -2, 2 } );
    k.push_back( std::vector< int >{ 1, 0, 0, 2, 0 } );
    k.push_back( std::vector< int >{ 1, 0, 2, -2, 1 } ); // 45
    k.push_back( std::vector< int >{ 0, 0, 0, -2, 1 } );
    k.push_back( std::vector< int >{ 0, -1, 2, -2, 1 } );
    k.push_back( std::vector< int >{ 2, 0, 2, 0, 1 } );
    k.push_back( std::vector< int >{ 1, -1, 0, 0, 0 } );
    k.push_back( std::vector< int >{ 1, 0, 0, -1, 0 } ); //50
    k.push_back( std::vector< int >{ 0, 0, 0, 1, 0 } );
    k.push_back( std::vector< int >{ 0, 1, 0, -2, 0 } );
    k.push_back( std::vector< int >{ 1, 0, -2, 0, 0 } );
    k.push_back( std::vector< int >{ 2, 0,  0, -2, 1 } );
    k.push_back( std::vector< int >{ 0, 1, 2, -2, 1 } ); // 55
    k.push_back( std::vector< int >{ 1, 1, 0, 0, 0 } );
    k.push_back( std::vector< int >{ 1, -1, 0, -1, 0 } );
    k.push_back( std::vector< int >{ -1, -1, 2, 2, 2 } );
    k.push_back( std::vector< int >{ 0, -1, 2, 2, 2 } );
    k.push_back( std::vector< int >{ 1, -1, 2, 0, 2 } ); //60
    k.push_back( std::vector< int >{ 3, 0, 2, 0, 2 } );
    k.push_back( std::vector< int >{ -2, 0, 2, 0, 2 } );
    k.push_back( std::vector< int >{ 1, 0, 2, 0, 0 } );
    k.push_back( std::vector< int >{ -1, 0, 2, 4, 2 } );
    k.push_back( std::vector< int >{ 1, 0, 0, 0, 2 } ); //65
    k.push_back( std::vector< int >{ -1, 0, 2, -2, 1 } );
    k.push_back( std::vector< int >{ 0, -2, 2, -2, 1 } );
    k.push_back( std::vector< int >{ -2, 0, 0, 0, 1 } );
    k.push_back( std::vector< int >{ 2, 0, 0, 0, 1 } );
    k.push_back( std::vector< int >{ 3, 0, 0, 0, 0 } ); // 70
    k.push_back( std::vector< int >{ 1, 1, 2, 0, 2 } );
    k.push_back( std::vector< int >{ 0, 0, 2, 1, 2 } );
    k.push_back( std::vector< int >{ 1, 0, 0, 2, 1 } );
    k.push_back( std::vector< int >{ 1, 0, 2, 2, 1 } );
    k.push_back( std::vector< int >{ 1, 1, 0, -2, 1 } ); // 75
    k.push_back( std::vector< int >{ 0, 1, 0, 2, 0 } );
    k.push_back( std::vector< int >{ 0, 1, 2, -2, 0 } );
    k.push_back( std::vector< int >{ 0, 1, -2, 2, 0 } );
    k.push_back( std::vector< int >{ 1, 0, -2, 2, 0 } );
    k.push_back( std::vector< int >{ 1, 0, -2, -2, 0 } ); //80
    k.push_back( std::vector< int >{ 1, 0, 2, -2, 0 } );
    k.push_back( std::vector< int >{ 1, 0, 0, -4, 0 } );
    k.push_back( std::vector< int >{ 2, 0, 0, -4, 0 } );
    k.push_back( std::vector< int >{ 0, 0, 2, 4, 2 } );
    k.push_back( std::vector< int >{ 0, 0, 2, -1, 2 } ); // 85
    k.push_back( std::vector< int >{ -2, 0, 2, 4, 2 } );
    k.push_back( std::vector< int >{ 2, 0, 2, 2, 2 } );
    k.push_back( std::vector< int >{ 0, -1, 2, 0, 1 } );
    k.push_back( std::vector< int >{ 0, 0, -2, 0, 1 } );
    k.push_back( std::vector< int >{ 0, 0, 4, -2, 2 } ); // 90
    k.push_back( std::vector< int >{ 0, 1, 0, 0, 2 } );
    k.push_back( std::vector< int >{ 1, 1, 2, -2, 2 } );
    k.push_back( std::vector< int >{ 3, 0, 2, -2, 2 } );
    k.push_back( std::vector< int >{ -2, 0, 2, 2, 2 } );
    k.push_back( std::vector< int >{ -1, 0, 0, 0, 2 } ); // 95
    k.push_back( std::vector< int >{ 0, 0, -2, 2, 1 } );
    k.push_back( std::vector< int >{ 0, 1, 2, 0, 1 } );
    k.push_back( std::vector< int >{ -1, 0, 4, 0, 2 } );
    k.push_back( std::vector< int >{ 2, 1, 0, -2, 0 } );
    k.push_back( std::vector< int >{ 2, 0, 0, 2, 0 } ); // 100
    k.push_back( std::vector< int >{ 2, 0, 2, -2, 1 } );
    k.push_back( std::vector< int >{ 2, 0, -2, 0, 1 } );
    k.push_back( std::vector< int >{ 1, -1, 0, -2, 0 } );
    k.push_back( std::vector< int >{ -1, 0, 0, 1, 1 } );
    k.push_back( std::vector< int >{ -1, -1, 0, 2, 1 } ); // 105
    k.push_back( std::vector< int >{ 0, 1, 0, 1, 0 } );

    A.push_back( std::vector< double >{ -171996, -174.2 } );
    A.push_back( std::vector< double >{ -13187, -1.6} );
    A.push_back( std::vector< double >{ -2274, -0.2} );
    A.push_back( std::vector< double >{ 2062, 0.2} );
    A.push_back( std::vector< double >{ -1426, 3.4} );
    A.push_back( std::vector< double >{ 712, 0.1} );
    A.push_back( std::vector< double >{ -517, 1.2} );
    A.push_back( std::vector< double >{ -386, -0.4 } );
    A.push_back( std::vector< double >{ -301, 0} );
    A.push_back( std::vector< double >{ 217, -0.5} ); // 10
    A.push_back( std::vector< double >{ 158, 0} );
    A.push_back( std::vector< double >{ 129, 0.1} );
    A.push_back( std::vector< double >{ 123, 0} );
    A.push_back( std::vector< double >{ 63, 0.1} );
    A.push_back( std::vector< double >{ 63, 0} );
    A.push_back( std::vector< double >{ -59, 0} );
    A.push_back( std::vector< double >{ -58, -0.1} );
    A.push_back( std::vector< double >{ -51, 0} );
    A.push_back( std::vector< double >{ -48, 0} );
    A.push_back( std::vector< double >{ 46, 0} );
    A.push_back( std::vector< double >{ -38, 0} );
    A.push_back( std::vector< double >{ -31, 0} );
    A.push_back( std::vector< double >{ 29, 0} );
    A.push_back( std::vector< double >{ 29, 0} );
    A.push_back( std::vector< double >{ 26, 0} );
    A.push_back( std::vector< double >{ -22, 0} );
    A.push_back( std::vector< double >{ 21, 0} );
    A.push_back( std::vector< double >{ 17, -0.1} );
    A.push_back( std::vector< double >{ -16, 0.1} );
    A.push_back( std::vector< double >{ 16, 0} );
    A.push_back( std::vector< double >{ -15, 0} );
    A.push_back( std::vector< double >{ -13, 0} );
    A.push_back( std::vector< double >{ -12, 0} );
    A.push_back( std::vector< double >{ 11, 0} );
    A.push_back( std::vector< double >{ -10, 0} );
    A.push_back( std::vector< double >{ -8, 0} );
    A.push_back( std::vector< double >{ -7, 0} );
    A.push_back( std::vector< double >{ -7, 0} );
    A.push_back( std::vector< double >{ -7, 0} );
    A.push_back( std::vector< double >{ 7, 0} ); // 40 
    A.push_back( std::vector< double >{ -6, 0 } );
    A.push_back( std::vector< double >{ -6, 0 } );
    A.push_back( std::vector< double >{ 6, 0 } );
    A.push_back( std::vector< double >{ 6, 0 } );
    A.push_back( std::vector< double >{ 6, 0 } ); // 45
    A.push_back( std::vector< double >{ -5, 0 } );
    A.push_back( std::vector< double >{ -5, 0 } );
    A.push_back( std::vector< double >{ -5, 0 } );
    A.push_back( std::vector< double >{ 5, 0 } );
    A.push_back( std::vector< double >{ -4, 0 } ); // 50
    A.push_back( std::vector< double >{ -4, 0 } );
    A.push_back( std::vector< double >{ -4, 0 } );
    A.push_back( std::vector< double >{ 4, 0 } );
    A.push_back( std::vector< double >{ 4, 0 } );
    A.push_back( std::vector< double >{ 4, 0 } ); // 55
    A.push_back( std::vector< double >{ -3, 0 } );
    A.push_back( std::vector< double >{ -3, 0 } );
    A.push_back( std::vector< double >{ -3, 0 } );
    A.push_back( std::vector< double >{ -3, 0 } );
    A.push_back( std::vector< double >{ -3, 0 } ); // 60
    A.push_back( std::vector< double >{ -3, 0 } );
    A.push_back( std::vector< double >{ -3, 0 } );
    A.push_back( std::vector< double >{ 3, 0 } );
    A.push_back( std::vector< double >{ -2, 0 } );
    A.push_back( std::vector< double >{ -2, 0 } ); // 65
    A.push_back( std::vector< double >{ -2, 0 } );
    A.push_back( std::vector< double >{ -2, 0 } );
    A.push_back( std::vector< double >{ -2, 0 } );
    A.push_back( std::vector< double >{ 2, 0 } );
    A.push_back( std::vector< double >{ 2, 0 } ); // 70
    A.push_back( std::vector< double >{ 2, 0 } );
    A.push_back( std::vector< double >{ 2, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } ); // 75
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } ); // 80
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } ); // 85
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ -1, 0 } );
    A.push_back( std::vector< double >{ 1, 0 } ); // 90
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); // 95 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); // 100
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); 
    A.push_back( std::vector< double >{ 1, 0 } ); // 105
    A.push_back( std::vector< double >{ 1, 0 } ); 

    B.push_back( std::vector< double >{ 92025, 8.9} );
    B.push_back( std::vector< double >{ 5736, -3.1} );
    B.push_back( std::vector< double >{ 977, -0.5} );
    B.push_back( std::vector< double >{ -895, 0.5} );
    B.push_back( std::vector< double >{ 54, -0.1} );
    B.push_back( std::vector< double >{ -7, 0} );
    B.push_back( std::vector< double >{ 224, -0.6} );
    B.push_back( std::vector< double >{ 200, 0} );
    B.push_back( std::vector< double >{ 129, -0.1} );
    B.push_back( std::vector< double >{ -95, 0.3} );
    B.push_back( std::vector< double >{ -1, 0} );
    B.push_back( std::vector< double >{ -70, 0} );
    B.push_back( std::vector< double >{ -53, 0} );
    B.push_back( std::vector< double >{ -33, 0} );
    B.push_back( std::vector< double >{ -2, 0} );
    B.push_back( std::vector< double >{ 26, 0} );
    B.push_back( std::vector< double >{ 32, 0} );
    B.push_back( std::vector< double >{ 27, 0} );
    B.push_back( std::vector< double >{ 1, 0} );
    B.push_back( std::vector< double >{ -24, 0} );
    B.push_back( std::vector< double >{ 16, 0 } );
    B.push_back( std::vector< double >{ 13, 0} );
    B.push_back( std::vector< double >{ -1, 0} );
    B.push_back( std::vector< double >{ -12, 0} );
    B.push_back( std::vector< double >{ -1, 0} );
    B.push_back( std::vector< double >{ 0, 0} );
    B.push_back( std::vector< double >{ -10, 0} );
    B.push_back( std::vector< double >{ 0, 0} );
    B.push_back( std::vector< double >{ 7, 0} );
    B.push_back( std::vector< double >{ -8, 0} );
    B.push_back( std::vector< double >{ 9, 0} );
    B.push_back( std::vector< double >{ 7, 0} );
    B.push_back( std::vector< double >{ 6, 0} );
    B.push_back( std::vector< double >{ 0, 0} );
    B.push_back( std::vector< double >{ 5, 0} );
    B.push_back( std::vector< double >{ 3, 0} );
    B.push_back( std::vector< double >{ 3, 0} );
    B.push_back( std::vector< double >{ 3, 0} );
    B.push_back( std::vector< double >{ 0, 0} );
    B.push_back( std::vector< double >{ -3, 0} ); // 40
    B.push_back( std::vector< double >{ 3, 0 } ); 
    B.push_back( std::vector< double >{ 3, 0 } ); 
    B.push_back( std::vector< double >{ -3, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ -3, 0 } ); //45
    B.push_back( std::vector< double >{ 3, 0 } ); 
    B.push_back( std::vector< double >{ 3, 0 } ); 
    B.push_back( std::vector< double >{ 3, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 50 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ -2, 0 } ); 
    B.push_back( std::vector< double >{ -2, 0 } ); // 55 
    B.push_back( std::vector< double >{ 0, 0  } ); 
    B.push_back( std::vector< double >{ 0, 0  } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); // 60 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); // 65 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ -1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 70
    B.push_back( std::vector< double >{ -1, 0 } ); 
    B.push_back( std::vector< double >{ -1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 75 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 80
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 85
    B.push_back( std::vector< double >{ 1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 90
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ -1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ -1, 0 } ); 
    B.push_back( std::vector< double >{ -1, 0 } ); // 95 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 100
    B.push_back( std::vector< double >{ -1, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); 
    B.push_back( std::vector< double >{ 0, 0 } ); // 105
    B.push_back( std::vector< double >{ 0, 0 } ); 
/*
    double c = pow(10,-4); 
    for( int i=0; i<A.size(); ++i )
    {
        for( int j=0; j<2; ++j )
        {
            A.at(i).at(j) = A.at(i).at(j)*c;
            B.at(i).at(j) = A.at(i).at(j)*c;
        }
    }*/
}

#endif // DEF_H
