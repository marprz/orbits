#ifndef DEF_H
#define DEF_H

#include <list>
#include <deque>
#include <array>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

typedef std::array< std::array< double, 3 >, 3 > Matrix33;
typedef std::vector< std::vector< double > > MatrixD;
typedef std::vector< std::vector< int > > MatrixI;
typedef std::vector< double > Vector;

double kmuE = 398600.4418;
double kRE = 6738;
double kRE2 = kRE*kRE;
double kJ2 = 0.00108263;
double kc = 299792.458; // speed of light [km/s]

void saveVectorToFile( const std::vector< std::vector< Vector > >& positions, const std::string& outName, const int& rows = 0 )
{
    std::cout << "Writing to file " << outName << "..." << std::endl;
    int satellites_nb = positions.size();
    int positions_nb;
    if( rows != 0 )
        positions_nb = rows;
    else 
        positions_nb = positions.at(0).size();
    std::cout << "satellites_nb = " << satellites_nb << ", positions_nb = " << positions_nb << std::endl;
    std::ofstream outFile;
    outFile.open( outName.c_str(), std::ios::out );
    outFile << "function p = " << outName << "\n";
    outFile << "p=zeros(" << satellites_nb << ", " << positions_nb << ", 3 );\n";
    for( int i=0; i<satellites_nb; ++i )
    {
        int temp_pos_nb = positions.at(i).size();
        if( temp_pos_nb > 0 )
        {
            std::cout << "i=" << i << ", temp_pos_nb=" << temp_pos_nb << std::endl;
            outFile << "p(" << i+1 << ",:,:) = [ ";
            for( int j=0; j<temp_pos_nb-1; ++j )
            {
                double x = ((positions.at(i)).at(j)).at(0);
                double y = ((positions.at(i)).at(j)).at(1);
                double z = ((positions.at(i)).at(j)).at(2);
                outFile << x << " " << y << " " << z << ";\n";
//                std::cout << x << " " << y << " " << z << ";\n";
             //   std::cout << "acceleration: " << norm(Acceleration( (positions.at(i)).at(j))) << std::endl;
            }
            outFile << ((positions.at(i)).at(temp_pos_nb-1)).at(0) << " " << ((positions.at(i)).at(temp_pos_nb-1)).at(1) << " " << ((positions.at(i)).at(temp_pos_nb-1)).at(2) << "];\n";
        }
    }
    outFile.close();
}

void saveVectorToFile( const std::vector< Vector >& positions, const std::string& outName, const double& h=0, const int& rows = 0 )
{
    std::cout << "saveVectorToFile: h=" << h << ", " << outName << std::endl;
    std::string fileName( "outputs/");
    fileName += outName;
    fileName += ".m";
    std::cout << "positions.size() = " << positions.size() << std::endl;
    int positions_nb;
    if( rows != 0 )
        positions_nb = rows;
    else 
        positions_nb = positions.size();
    std::fstream outFile;
    outFile.open( fileName.c_str(), std::ios::out );
    outFile << "function p = " << outName << "\n";
    int col = 3;
    if( h != 0 )
        ++col;
    outFile << "p=zeros(" << positions_nb << ", " << col << " );\n";
    int temp_pos_nb = positions.size();
    if( temp_pos_nb > 0 )
    {
        outFile << "p = [ ";
        for( int j=0; j<temp_pos_nb; ++j )
       // for( int j=0; j<temp_pos_nb-1; ++j )
        {
            double x = (positions.at(j)).at(0);
            double y = (positions.at(j)).at(1);
            double z = (positions.at(j)).at(2);
            if( h == 0 )
            {
                if( temp_pos_nb-1 == j )
                    outFile << x << " " << y << " " << z << "];\n";
                else
                    outFile << x << " " << y << " " << z << ";\n";
            }
            else
            {
                double t = (j-4)*h;
                if( temp_pos_nb-1 == j )
                    outFile << t << " " << x << " " << y << " " << z << "];\n";
                else
                    outFile << t << " " << x << " " << y << " " << z << ";\n";
            }
 
        }
//        outFile << (positions.at(temp_pos_nb-1)).at(0) << " " << (positions.at(temp_pos_nb-1)).at(1) << " " << (positions.at(temp_pos_nb-1)).at(2) << "];\n";

    }
    outFile.close();
}



MatrixD createMatrix( std::vector< double > v )
{
    MatrixD M; 
    Vector v1;
    v1.push_back( v.at(0) );
    v1.push_back( v.at(1) );
    v1.push_back( v.at(2) );
    Vector v2;
    v2.push_back( v.at(3) );
    v2.push_back( v.at(4) );
    v2.push_back( v.at(5) );
    Vector v3;
    v3.push_back( v.at(6) );
    v3.push_back( v.at(7) );
    v3.push_back( v.at(8) );
    M.push_back( v1 );
    M.push_back( v2 );
    M.push_back( v3 );

    return M;
}

double norm( const Vector& v )
{
    double ret = 0;
    for( int i=0; i<v.size(); ++i )
        ret += (v.at(i))*(v.at(i));
    ret = sqrt( ret );
    return ret;
}

double VelocityMagnitude( const Vector& v )
{
    double mag = norm( v );
    double ra = kRE+20394; // apogee [km]
    return sqrt( kmuE*( 2/mag - 1/ra ) );
}

void print( Vector V, std::string name = "Vector" )
{
    if( name != "" )
    {
        std::cout << name << ": ";
    }
//    double r = sqrt( pow(V.at(0),2) + pow(V.at(1),2) + pow(V.at(2),2) );
    for( int i=0; i<V.size(); ++i )
    {
        std::setprecision(10);
        std::cout.width( 10 );
        std::cout << V.at(i) << " ";
    }
    std::cout << "   [ " << norm(V) << " ] " << std::endl;
}


void print( std::vector< Vector > v, std::string name = "Vector" )
{
    std::cout << name << ": " << std::endl;
    std::vector< Vector >::iterator it;
    for( it = v.begin(); it != v.end(); ++it )
        print( *it, "" );
    std::cout << std::endl;
}

void print( const std::deque< Vector >& v, std::string name = "Vector" )
{
    std::cout << name << ": " << std::endl;
    std::deque< Vector >::const_iterator it;
    for( it = v.begin(); it != v.end(); ++it )
        print( *it, "" );
    std::cout << std::endl;
}


Vector operator*( double d, Vector v )
{
    for( int i=0; i<v.size(); ++i )
        v.at(i) *= d;
    return v;
}

Vector operator*( Vector v, double d )
{
    return d*v;
}

Vector operator-( Vector v1, Vector v2 )
{
    if( v1.size() != v2.size() )
        std::cout << "ERROR!Vector operator-( Vector v1, Vector v2 ) " << std::endl;
    for( int i=0; i<v1.size(); ++i )
        v1.at(i) -= v2.at(i);
    return v1;
}

std::vector< Vector > operator-( std::vector< Vector > v )
{
    std::vector< Vector >::iterator it;
    for( it = v.begin(); it != v.end(); ++it )
        (*it ) = (-1)*(*it);
    return v;
}

Vector operator-( Vector v )
{
    return (-1)*v;
}

std::vector< Vector > operator-( std::vector< Vector > v1, const std::vector< Vector >& v2 )
{
    if( v1.size() != v2.size() )
        std::cout << "ERROR! Cannot calculate vector< Vector >-vector< Vector > => different sizes!" << std::endl;
    for( int i=0; i<v1.size(); ++i )
        v1.at(i) = v1.at(i)-v2.at(i);
    return v1;
}

Vector operator/( Vector v, double d )
{
    for( int i=0; i<v.size(); ++i )
        v.at(i) /= d;
    return v;
}

// mnożenie dwóch wektorów. 
// jeśli mają różne długości - mnożenie ostatnich elementów
double operator*( Vector v1, Vector v2 )
{
    double sum = 0;
    int s1 = v1.size();
    int s2 = v2.size();
    int m = s1 < s2 ? s1 : s2; 
    for( int i=0; i<m; ++i )
        sum += v1.at(s1-1-i)*v2.at(s2-1-i);
    return sum;
}

Vector operator+( double d, Vector v )
{
    for( int i=0; i<v.size(); ++i )
        v.at(i) += d;
    return v;
}

Vector operator+( const Vector& v1, const Vector& v2 )
{
    int m = v1.size() < v2.size() ? v1.size() : v2.size();
    Vector ret = v1.size() > v2.size() ? v1 : v2;

    for( int i=0; i<m; ++i )
        ret.at(i) = v1.at(i) + v2.at(i);

    return ret;
}

void print( MatrixD, std::string );
/*
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
    Vector v1 = { A/detM, D/detM, G/detM };
    Vector v2 = { B/detM, E/detM, H/detM };
    Vector v3 = { C/detM, F/detM, I/detM };

    MatrixD inverse = { v1, v2, v3 };
    return inverse;
}*/

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

/*
// TODO przerobic drukowanie na template
void print( MatrixD M, std::string name = "Matrix" )
{
    std::cout << std::endl << name << ": " << std::endl;
    std::cout << M.at(0).at(0) << " " << M.at(0).at(1) << " " << M.at(0).at(2) << std::endl;
    std::cout << M.at(1).at(0) << " " << M.at(1).at(1) << " " << M.at(1).at(2) << std::endl;
    std::cout << M.at(2).at(0) << " " << M.at(2).at(1) << " " << M.at(2).at(2) << std::endl;
}*/

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

    Vector v1;
    v1.push_back( 1 );
    v1.push_back( 0 );
    v1.push_back( 0 );
   
    Vector v2;
    v2.push_back( 0 );
    v2.push_back( cos(m) );
    v2.push_back( sin(m) );

    Vector v3;
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

    Vector v1;
    v1.push_back( cos(m) );
    v1.push_back( 0 );
    v1.push_back( -sin(m) );
   
    Vector v2;
    v2.push_back( 0 );
    v2.push_back( 1 );
    v2.push_back( 0 );

    Vector v3;
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

    Vector v1;
    v1.push_back( cos(m) );
    v1.push_back( sin(m) );
    v1.push_back( 0 );
   
    Vector v2;
    v2.push_back( -sin(m) );
    v2.push_back( cos(m) );
    v2.push_back( 0 );

    Vector v3;
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
        Vector v;
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

Vector multiply( MatrixD M, Vector v )
{
    Vector ret;
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

/*
Vector operator- ( const Vector& v1, const Vector& v2 )
{
    Vector ret;
    for( int i=0; i<3; ++i )
        ret.push_back( v1.at(i)-v2.at(i) );
    return ret;
}*/

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

void HairerParameters( Vector& ALPHA, std::vector< Vector >& A, Vector& AL )
{
    ALPHA.resize(11);
    ALPHA.at(0) = 0;
    ALPHA.at(1) = 8.48093651182655889791*0.01;
    ALPHA.at(2) = 1.69618730236531177958*0.1;
    ALPHA.at(3) = 8.82527661964732346426*0.1;
    ALPHA.at(4) = 3.57384241759677451843*0.1;
    ALPHA.at(5) = 1.17472338035267653574*0.1;
    ALPHA.at(6) = 6.42615758240322548157*0.1;
    ALPHA.at(7) = ALPHA.at(3);
    ALPHA.at(8) = ALPHA.at(2);
    ALPHA.at(9) = 0;
    ALPHA.at(10) = 1;

    A.resize(12);
    for( int i=0; i<12; ++i )
    {
        A.at(i).resize(11);
        for( int j=0; j<11; ++j )
            A.at(i).at(j) = 0;
    }
    
    A.at(1).at(0)  = 3.59631420588164200966*0.001;
    A.at(2).at(0)  = 4.79508560784218934622*0.001;
    A.at(2).at(1)  = 9.59017121568437869244*0.001;
    A.at(3).at(0)  = 1.87733546050627049772;
    A.at(3).at(1)  =-4.32661231970111510045;
    A.at(3).at(2)  = 2.83870439626131304487;
    A.at(4).at(0)  = 1.94709467473823596644*0.01;
    A.at(4).at(1)  = 0;
    A.at(4).at(2)  = 4.42810873502123510231*0.01;
    A.at(4).at(3)  = 1.09714031475079632769*0.0001;
    A.at(5).at(0)  = 4.64206697543819328333*0.001;
    A.at(5).at(1)  = 0;
    A.at(5).at(2)  = 2.93579155497804175774*0.001;
    A.at(5).at(3)  = 2.76435808114904080285*0.00001;
    A.at(5).at(4)  =-7.05627009491629727298*0.0001;
    A.at(6).at(0)  =-6.96029489417689828096*0.001;
    A.at(6).at(1)  = 0;
    A.at(6).at(2)  =-1.90054512658393594723*0.1;
    A.at(6).at(3)  = 7.68088148461498565180*0.1;
    A.at(6).at(4)  = 1.18704223034434343286*0.1;
    A.at(6).at(5)  = 2.84020002739066989630*0.1;
    A.at(7).at(0)  = 8.26380420350248571364*0.01;
    A.at(7).at(1)  = 0;
    A.at(7).at(2)  = 3.60887996724131292687*0.1;
    A.at(7).at(3)  = 1.01208078941331907108*0.001;
    A.at(7).at(4)  = 6.79462118331253280025*0.01;
    A.at(7).at(5)  =-2.04227957300564960765*0.1;
    A.at(7).at(6)  = 8.11711629853386060148*0.01;
    A.at(8).at(0)  = 2.38468700370378669757*0.01;
    A.at(8).at(1)  = 9.59017121568437869244*0.001;
    A.at(8).at(2)  = 1.52655867995367678583*0.1;
    A.at(8).at(3)  =-5.41908820084510828926*0.001;
    A.at(8).at(4)  =-3.93767256187270659404*0.01;
    A.at(8).at(5)  =-1.40338104799298795297*0.1; 
    A.at(8).at(6)  = 1.00094410252962054139*0.01;
    A.at(8).at(7)  = 3.41682516901140790025*0.001;
    A.at(9).at(0) =-4.25990467506340701948*0.01;
    A.at(9).at(1) = 0;
    A.at(9).at(2) =-5.93485416695720137147*0.1;
    A.at(9).at(3) = 9.81481906263695126914*0.001;
    A.at(9).at(4) = 7.13421080958931560696*0.01;
    A.at(9).at(5) = 2.91368066505515770072*0.1;
    A.at(9).at(6) =-1.57502473657776550315*0.01;
    A.at(9).at(7) =-6.97794897116397988241*0.001;
    A.at(9).at(8) = 2.86287666119249964845*0.1;
    A.at(10).at(0) =-6.14180937931931543943*0.1;
    A.at(10).at(1) = 0;
    A.at(10).at(2) =-1.75100055986614758156;
    A.at(10).at(3) =-6.76603738760857112582*0.001;
    A.at(10).at(4) = 3.95581030796451903257*0.1;
    A.at(10).at(5) = 1.13603914070927980965;
    A.at(10).at(6) = 5.64165276602374526099*0.01;
    A.at(10).at(7) = 3.58989223779823259146*0.01;
    A.at(10).at(8) = 7.48011913641736205396*0.1;
    A.at(10).at(9)= 5*0.1;
    A.at(11).at(0) = 0;
    A.at(11).at(1) = 0;
    A.at(11).at(2) =-8.30381269763468822042*0.01;
    A.at(11).at(3) = 0;
    A.at(11).at(4) = 1.78280368337326917339*0.1;
    A.at(11).at(5) = 1.67007309146871573764*0.1;
    A.at(11).at(6) = 9.91488201804162591694*0.01;
    A.at(11).at(7) = 2.22301690020519163947*0.01;
    A.at(11).at(8) = 8.30381269763468822042*0.01;
    A.at(11).at(9)= 3.33333333333333333333*0.01;
    A.at(11).at(10)= 0;
    
    AL.resize(11);
    AL.at(0) = 0;
    AL.at(1) = 0;
    AL.at(2) =-1*0.1;
    AL.at(3) = 0;
    AL.at(4) = 2.77429188517743176508*0.1;
    AL.at(5) = 1.89237478148923490158*0.1;
    AL.at(6) = 2.77429188517743176508*0.1;
    AL.at(7) = 1.89237478148923490158*0.1;
    AL.at(8) = 1*0.1;
    AL.at(9)= 3.33333333333333333333*0.01;
    AL.at(10)= 3.33333333333333333333*0.01;
}

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
