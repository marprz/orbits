#include <algorithm>
#include <numeric>
#include "gtest/gtest.h"
#include "Taylor.cpp"

TEST( TestDef, createMatrixTest )
{
    std::vector< double > v;
    v.push_back( 1 );
    v.push_back( 2 );
    v.push_back( 3 );
    v.push_back( 4 );
    v.push_back( 5 );
    v.push_back( 6 );
    v.push_back( 7 );
    v.push_back( 8 );
    v.push_back( 9 );
    MatrixD out = createMatrix( v );

    ASSERT_EQ( out.at(0).at(0), 1 );
    ASSERT_EQ( out.at(0).at(1), 2 );
    ASSERT_EQ( out.at(0).at(2), 3 );
    ASSERT_EQ( out.at(1).at(0), 4 );
    ASSERT_EQ( out.at(1).at(1), 5 );
    ASSERT_EQ( out.at(1).at(2), 6 );
    ASSERT_EQ( out.at(2).at(0), 7 );
    ASSERT_EQ( out.at(2).at(1), 8 );
    ASSERT_EQ( out.at(2).at(2), 9 );
}

TEST( TestDef, TranspositionTest )
{
    MatrixD M = { { 1, 2, 3}, { 4, 5, 6 }, { 7, 8, 9 } };
    Transpose( M );
    ASSERT_EQ( M.at(0).at(0), 1 );
    ASSERT_EQ( M.at(0).at(1), 4 );
    ASSERT_EQ( M.at(0).at(2), 7 );
    ASSERT_EQ( M.at(1).at(0), 2 );
    ASSERT_EQ( M.at(1).at(1), 5 );
    ASSERT_EQ( M.at(1).at(2), 8 );
    ASSERT_EQ( M.at(2).at(0), 3 );
    ASSERT_EQ( M.at(2).at(1), 6 );
    ASSERT_EQ( M.at(2).at(2), 9 );

}

TEST( TestDef, R1Test )
{
    double angle = 1.570796325; // [rad] = 90 degrees
    MatrixD r1 = R1( angle );
    
    double delta = 0.0001;
    ASSERT_EQ( r1.at(0).at(0), 1 );
    ASSERT_EQ( r1.at(0).at(1), 0 );
    ASSERT_EQ( r1.at(0).at(2), 0 );
    ASSERT_EQ( r1.at(1).at(0), 0 );
    ASSERT_GE( r1.at(1).at(1), -delta );
    ASSERT_LE( r1.at(1).at(1), delta );
    ASSERT_GE( r1.at(1).at(2), 1-delta );
    ASSERT_LE( r1.at(1).at(2), 1+delta );
    ASSERT_EQ( r1.at(2).at(0), 0 );
    ASSERT_GE( r1.at(2).at(1), -1-delta );
    ASSERT_LE( r1.at(2).at(1), -1+delta );
    ASSERT_GE( r1.at(2).at(2), -delta );
    ASSERT_LE( r1.at(2).at(2), delta );
}

TEST( TestDef, R2Test )
{
    double OneDegToRad = 0.0174532925;
    double angle = 90*OneDegToRad; // [rad] = 90 degrees
    MatrixD r2 = R2( angle );
    double delta = 0.0001;

    ASSERT_GE( r2.at(0).at(0), -delta );
    ASSERT_LE( r2.at(0).at(0), delta );
    ASSERT_EQ( r2.at(0).at(1), 0 );
    ASSERT_GE( r2.at(0).at(2), -1-delta );
    ASSERT_LE( r2.at(0).at(2), -1+delta );
    ASSERT_EQ( r2.at(1).at(0), 0 );
    ASSERT_EQ( r2.at(1).at(1), 1 );
    ASSERT_EQ( r2.at(1).at(2), 0 );
    ASSERT_GE( r2.at(2).at(0), 1-delta );
    ASSERT_LE( r2.at(2).at(0), 1+delta );
    ASSERT_EQ( r2.at(2).at(1), 0 );
    ASSERT_GE( r2.at(2).at(2), -delta );
    ASSERT_LE( r2.at(2).at(2), delta );
}

TEST( TestDef, R3Test )
{
    double OneDegToRad = 0.0174532925;
    double angle = 90*OneDegToRad; // [rad] = 90 degrees
    MatrixD r3 = R3( angle );
    double delta = 0.0001;

    ASSERT_GE( r3.at(0).at(0), -delta );
    ASSERT_LE( r3.at(0).at(0), delta );
    ASSERT_GE( r3.at(0).at(1), 1-delta );
    ASSERT_LE( r3.at(0).at(1), 1+delta );
    ASSERT_EQ( r3.at(0).at(2), 0 );
    ASSERT_GE( r3.at(1).at(0), -1-delta );
    ASSERT_LE( r3.at(1).at(0), -1+delta );
    ASSERT_GE( r3.at(1).at(1), -delta );
    ASSERT_LE( r3.at(1).at(1), delta );
    ASSERT_EQ( r3.at(1).at(2), 0 );
    ASSERT_EQ( r3.at(2).at(0), 0 );
    ASSERT_EQ( r3.at(2).at(1), 0 );
    ASSERT_EQ( r3.at(2).at(2), 1 );
}

TEST( TestDef, multiplicationTest1 )
{
    std::vector< double > v1 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    std::vector< double > v2 = { 2, 0, 0, 0, 2, 0, 0, 0, 2 };
    MatrixD mIn1 = createMatrix( v1 );
    MatrixD mIn2 = createMatrix( v2 );
    MatrixD mOut = multiply( mIn1, mIn2 );
    ASSERT_EQ( mOut.at(0).at(0), 2 );
    ASSERT_EQ( mOut.at(0).at(1), 0 );
    ASSERT_EQ( mOut.at(0).at(2), 0 );
    ASSERT_EQ( mOut.at(1).at(0), 0 );
    ASSERT_EQ( mOut.at(1).at(1), 2 );
    ASSERT_EQ( mOut.at(1).at(2), 0 );
    ASSERT_EQ( mOut.at(2).at(0), 0 );
    ASSERT_EQ( mOut.at(2).at(1), 0 );
    ASSERT_EQ( mOut.at(2).at(2), 2 );
}

TEST( TestDef, multiplicationTest2 )
{
    std::vector< double > v1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    std::vector< double > v2 = { 1, 1, 1, 2, 2, 2, 3, 3, 3 };
    MatrixD mIn1 = createMatrix( v1 );
    MatrixD mIn2 = createMatrix( v2 );
    MatrixD mOut = multiply( mIn1, mIn2 );
    ASSERT_EQ( mOut.at(0).at(0), 14 );
    ASSERT_EQ( mOut.at(0).at(1), 14 );
    ASSERT_EQ( mOut.at(0).at(2), 14 );
    ASSERT_EQ( mOut.at(1).at(0), 32 );
    ASSERT_EQ( mOut.at(1).at(1), 32 );
    ASSERT_EQ( mOut.at(1).at(2), 32 );
    ASSERT_EQ( mOut.at(2).at(0), 50 );
    ASSERT_EQ( mOut.at(2).at(1), 50 );
    ASSERT_EQ( mOut.at(2).at(2), 50 );
}

TEST( TestDef, multiplicationVectorTest1 )
{
    std::vector< double > v1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    std::vector< double > v2 = { 2, 3, 4 };
    MatrixD mIn1 = createMatrix( v1 );
    Vector vOut = multiply( mIn1, v2 );
    ASSERT_EQ( vOut.at(0), 20 );
    ASSERT_EQ( vOut.at(1), 47 );
    ASSERT_EQ( vOut.at(2), 74 );
}

TEST( TestDef, DaysInYearTest )
{
    Time t1( 2000, 1, 1, 0, 1, 0 );
    ASSERT_EQ( daysInYear( t1 ), 1 );
    Time t2( 1999, 12, 31, 0, 1, 0 );
    ASSERT_EQ( daysInYear( t2 ), 365 );
    Time t3( 1989, 3, 31, 10, 1, 0 );
    ASSERT_EQ( daysInYear( t3 ), 90 );
}

TEST( TestDef, LeapYearsSince1980 )
{
    ASSERT_EQ( leapYearsSince1980( 1980 ), 0 );
    ASSERT_EQ( leapYearsSince1980( 1984 ), 0 );
    ASSERT_EQ( leapYearsSince1980( 1985 ), 1 );
    ASSERT_EQ( leapYearsSince1980( 1988 ), 1 );
    ASSERT_EQ( leapYearsSince1980( 1993 ), 3 );
}

TEST( TestDef, JulianDateTest )
{
    ASSERT_EQ( JulianDate( 1950, 1, 1, 0 ), 2433282.5 );
    ASSERT_EQ( JulianDate( 2010, 1, 1, 0 ), 2455197.5 );
    Time t1( 2000, 1, 1, 12, 0, 0 );
    double expJD = 2451545;
    ASSERT_EQ( JulianDate( t1 ), expJD );
}

TEST( TestDef, ArcsecToRadTest )
{
    double r = 0.0000048481368111; 
    double eps = pow(10,-10);
    ASSERT_GE( arcsecToRad( 1 ), r-eps );
    ASSERT_LE( arcsecToRad( 1 ), r+eps );
}

TEST( TestDef, LatLonToXYZTest )
{
    double eps = 1; // 1 [m]
    std::vector< double > latlon1 = { 26500000, 0, 0 };
    std::vector< double > outXYZ1 = LatLonToXYZ( latlon1 );
    std::vector< double > expXYZ1 = { 26500000, 0, 0};
    Vector v1 = outXYZ1-expXYZ1;
    ASSERT_LE( norm(v1), eps );

    std::vector< double > latlon2 = { 26500000, 90, 0 };
    std::vector< double > outXYZ2 = LatLonToXYZ( latlon2 );
    std::vector< double > expXYZ2 = { 0, 0, 26500000};;
    Vector v2 = outXYZ2-expXYZ2;
    ASSERT_LE( norm(v2), eps );

    std::vector< double > latlon3 = { 26500000, 0, 180 };
    std::vector< double > outXYZ3 = LatLonToXYZ( latlon3 );
    std::vector< double > expXYZ3 = { -26500000, 0, 0};
    Vector v3 = outXYZ3-expXYZ3;
    ASSERT_LE( norm(v3), eps );

}

TEST( TestDef, JulianSince2000Century )
{
    // 1 == 100 years
    // 0.01 == 1 year
    double eps = 0.0001;
    Time t1( 2000, 1, 1, 12, 0, 0 );
    ASSERT_EQ( JulianSince2000Century( t1 ), 0 );
    Time t2( 2050, 1, 1, 12, 0, 0 );
    ASSERT_GE( JulianSince2000Century( t2 ), 0.5-eps );
    ASSERT_LE( JulianSince2000Century( t2 ), 0.5+eps );
    Time t3( 2100, 1, 1, 12, 0, 0 );
    ASSERT_EQ( JulianSince2000Century( t3 ), 1 );
    Time t4( 2010, 7, 1, 12, 0, 0 );
    ASSERT_GE( JulianSince2000Century( t4 ), 0.105-eps );
    ASSERT_LE( JulianSince2000Century( t4 ), 0.105+eps );
}

TEST( TestTaylor, factorial )
{
    ASSERT_EQ( fac(1), 1 );
    ASSERT_EQ( fac(2), 2 );
    ASSERT_EQ( fac(3), 6 );
    ASSERT_EQ( fac(4), 24 );
    ASSERT_EQ( fac(5), 120 );
    ASSERT_EQ( fac(6), 720 );
    ASSERT_EQ( fac(8), 40320 );
}

TEST( TestTaylor, serie1 )
{
    // f(x) = sin(x), x=0
    double h = 1;
    std::vector< double > v;
    v.push_back( 0 );
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );
    v.push_back( 0 );
    v.push_back( 1 );
    //EXPECT_EQ( taylorSerie( h, v ), 0.84147 ); 
    EXPECT_LE( taylorSerie( h, v ), 0.842 );  // less or equal
    EXPECT_GE( taylorSerie( h, v ), 0.841 ); // greater or equal
}

TEST( TestTaylor, serie2 )
{
    // f(x) = cos(x), x = pi
    double h = 3.141592654/3;
    double h2 = h*h;
    double h4 = h2*h2;
    double h6 = h4*h2;
    double delta = 0.001;
    std::vector< double > v;
    v.push_back( 1 );
    EXPECT_EQ( taylorSerie( h, v ), 1 );
    v.push_back( 0 );
    EXPECT_EQ( taylorSerie( h, v ), 1 );
    v.push_back( -1 );
    EXPECT_EQ( taylorSerie( h, v ), 1-(h2)/2 );
    v.push_back( 0 );
    EXPECT_EQ( taylorSerie( h, v ), 1-(h2)/2 );
    v.push_back( 1 );
    EXPECT_EQ( taylorSerie( h, v ), 1-(h2)/2+(h4)/24 );
    v.push_back( 0 );
    EXPECT_EQ( taylorSerie( h, v ), 1-(h2)/2+(h4)/24 );
    v.push_back( -1 );
    EXPECT_EQ( taylorSerie( h, v ), 1-(h2)/2+(h4)/24-(h6)/720 );
    EXPECT_LE( taylorSerie( h, v ), cos(h)+delta );
    EXPECT_GE( taylorSerie( h, v ), cos(h)-delta );
}

TEST( TestTaylor, series1 )
{
    // f(x) = cos(x), szereg: f(-2)+f(-1.5)+f(-1)+...+f(1.5)+f(2)
    int n = 9;
    double h = 0.5;
    std::vector< double > v;
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );
    v.push_back( 0 );
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );
    v.push_back( 0 );
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );

    std::vector< double > out = taylorSeries( n, h, v );
    double sum = 0;
    for(int i=0; i<out.size(); ++i )
    {
        sum += out.at(i);
    }

    EXPECT_GE( sum, 3.144 );
    EXPECT_LE( sum, 3.145 );
}

TEST( TestTaylor, series2 )
{
    // f(x) = sin(x), szereg: f(-2)+f(-1.5)+f(-1)+...+f(1.5)+f(2)
    int n = 9;
    double h = 0.5;
    std::vector< double > v;
    v.push_back( 0 );
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );
    v.push_back( 0 );
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );
    v.push_back( 0 );
    v.push_back( 1 );
    v.push_back( 0 );
    v.push_back( -1 );

    std::vector< double > out = taylorSeries( n, h, v );
    double sum = 0;
    for(int i=0; i<out.size(); ++i )
    {
        sum += out.at(i);
    }

    EXPECT_EQ( sum, 0 );
    EXPECT_GE( sum, -0.001 );
    EXPECT_LE( sum, 0.001 );
}
int main( int argc, char **argv )
{
    testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();
}
