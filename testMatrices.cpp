#include <algorithm>
#include <numeric>
#include <tuple>
#include "gtest/gtest.h"
#include "Conversion.cpp"
#include "Def.h"

TEST( TestConversion, MoonMeanAnomalyTest )
{
    Time t( 2010, 5, 1, 0, 0, 0);
    double T = JulianSince2000Century( t );
    double expJD = 2455317.5;
    ASSERT_EQ( JulianDate( t ), expJD );

    Time t2000( 2000, 1, 1, 12, 0, 0);
    double jd2000 = JulianDate( t ) - JulianDate( t2000 );

    Time t2( 2015, 3, 23, 0, 0, 0 );
    double eps = 0.000001;
    ASSERT_GE( JulianSince2000Century( t2 ), 0.152211-eps );
    ASSERT_LE( JulianSince2000Century( t2 ), 0.152211+eps );

    Time t3( 2011, 11, 11, 11, 11, 11 );
    ASSERT_GE( JulianSince2000Century( t3 ), 0.118603-eps );
    ASSERT_LE( JulianSince2000Century( t3 ), 0.118603+eps );

    //www.stjarnhimlen.se/comp/tutorial.html
    Time t4( 1990, 4, 19, 0, 0, 0 );
    double expAnomaly4 = -46173.9046;
    double eps4 = 0.01;
    ASSERT_GE( MoonMeanAnomaly( t4 ), expAnomaly4-eps4 );
    ASSERT_LE( MoonMeanAnomaly( t4 ), expAnomaly4+eps4 );
}

TEST( TestConversion, SunMeanAnomalyTest )
{
    //www.stjarnhimlen.se/comp/tutorial.html
    Time t( 1990, 4, 19, 0, 0, 0 );
    double T = JulianSince2000Century( t );
    double expAnomaly = -3135.9347;
/*    expAnomaly = fmod( expAnomaly, 360 );
    if( expAnomaly < 0 )
        expAnomaly += 360;*/
    double eps = 0.01;
    ASSERT_GE( SunMeanAnomaly( T ), expAnomaly-eps );
    ASSERT_LE( SunMeanAnomaly( T ), expAnomaly+eps );
}

TEST( TestConversion, MoonMeanArgOfLatitudeTest )
{
    Time t( 1990, 4, 19, 0, 0, 0 );
    double T = JulianSince2000Century( t );
    double exp_mean_latitude = 1.84; // deg
    double calc_mean_latitude = fmod( MoonMeanArgOfLatitude( T ), 360 );
    if( calc_mean_latitude < 0 )
        calc_mean_latitude += 360;
    double eps = 0.01;
//    std::cout << "MoonMeanArgOfLatitude: " << fmod(MoonMeanArgOfLatitude( T ),360)+360 << std::endl;
    ASSERT_GE( calc_mean_latitude, exp_mean_latitude-eps );
    ASSERT_LE( calc_mean_latitude, exp_mean_latitude+eps );
}

// D
TEST( TestConversion, MoonMeanElongationTest )
{   
    Time t( 1990, 4, 19, 0, 0, 0 );
    double T = JulianSince2000Century( t );
    double calc_elongation =  MoonMeanElongation( T );
    double exp_elongation = 287.7401;
    calc_elongation = Mod360( calc_elongation );
    double eps = 0.01;
    ASSERT_GE( calc_elongation, exp_elongation-eps );
    ASSERT_LE( calc_elongation, exp_elongation+eps );
}

// Omega
TEST( TestConversion, MeanLongitudeOfAscendingLunarNodeTest )
{
    Time t( 1990, 4, 19, 0, 0, 0 );
    double T = JulianSince2000Century( t );
    double calc_longitude = MeanLongitudeOfAscendingLunarNode( T );
    double exp_longitude = 312.7381;
    double eps = 0.01;
    ASSERT_GE( calc_longitude, exp_longitude-eps );
    ASSERT_LE( calc_longitude, exp_longitude+eps );

}

TEST( TestConversion, EpsilonTest1 )
{
    Time t( 1990, 4, 19, 0, 0, 0 );
    double exp_epsilon = 23+26/60+25.991/3600; 
    double eps = 0.5;
    ASSERT_GE( Epsilon( t ), exp_epsilon-eps );
    ASSERT_LE( Epsilon( t ), exp_epsilon+eps );
}

TEST( TestConversion, EpsilonTest2 )
{
    Time t( 2010, 5, 5, 0, 0, 0 );
    double exp_epsilon = 23+26/60+16.608/3600; 
    double eps = 0.5;
    ASSERT_GE( Epsilon( t ), exp_epsilon-eps );
    ASSERT_LE( Epsilon( t ), exp_epsilon+eps );
}

TEST( TestConversion, NutationParametersTest1 )
{
    Time t( 1990, 4, 19, 0, 0, 0 );
    double eps = 0.001; // deg
    auto tup = NutationParameters( t );
    double c_dpsi = std::get< 0 >( tup );
    double c_deps = std::get< 1 >( tup );
    std::array< double, 5 > c_alpha = std::get< 2 >( tup );

    double exp_deps = 6.561/3600; // deg
    ASSERT_GE( c_deps, exp_deps-eps );
    ASSERT_LE( c_deps, exp_deps+eps );

    double exp_dpsi = 11.632/3600; // deg
    ASSERT_GE( c_dpsi, exp_dpsi-eps );
    ASSERT_LE( c_dpsi, exp_dpsi+eps );
}

TEST( TestConversion, NutationParametersTest2 )
{
    Time t( 2010, 5, 5, 0, 0, 0 );
    double eps = 0.001; // deg
    auto tup = NutationParameters( t );
    double c_dpsi = std::get< 0 >( tup );
    double c_deps = std::get< 1 >( tup );
    std::array< double, 5 > c_alpha = std::get< 2 >( tup );

    double exp_deps = 2.456/3600; // deg
    ASSERT_GE( c_deps, exp_deps-eps );
    ASSERT_LE( c_deps, exp_deps+eps );

    double exp_dpsi = 15.579/3600; // deg
    ASSERT_GE( c_dpsi, exp_dpsi-eps );
    ASSERT_LE( c_dpsi, exp_dpsi+eps );
}
TEST( TestMatrices, PrecessionTest1 )
{
    Time t( 2010, 5, 5, 0, 0, 0 );
    MatrixD P = precession( t );
    double eps = 0.000001;
    ASSERT_GE( P.at(0).at(0), 0.9999968-eps );
    ASSERT_LE( P.at(0).at(0), 0.9999968+eps );
    ASSERT_GE( P.at(1).at(0), 0.0023120-eps );
    ASSERT_LE( P.at(1).at(0), 0.0023120+eps );
    ASSERT_GE( P.at(2).at(0), 0.0010046-eps );
    ASSERT_LE( P.at(2).at(0), 0.0010046+eps );
    ASSERT_GE( P.at(0).at(1),-0.0023120-eps );
    ASSERT_LE( P.at(0).at(1),-0.0023120+eps );
    ASSERT_GE( P.at(1).at(1), 0.9999973-eps );
    ASSERT_LE( P.at(1).at(1), 0.9999973+eps );
    ASSERT_GE( P.at(2).at(1),-0.0000012-eps );
    ASSERT_LE( P.at(2).at(1),-0.0000012+eps );
    ASSERT_GE( P.at(0).at(2),-0.0010045-eps );
    ASSERT_LE( P.at(0).at(2),-0.0010045+eps );
    ASSERT_GE( P.at(1).at(2),-0.0000011-eps );
    ASSERT_LE( P.at(1).at(2),-0.0000011+eps );
    ASSERT_GE( P.at(2).at(2), 0.9999995-eps );
    ASSERT_LE( P.at(2).at(2), 0.9999995+eps );
}

TEST( TestMatrices, NutationTest1 )
{
    Time t( 2010, 5, 5, 0, 0, 0 );
    MatrixD N = nutation( t );
    double eps = 0.0001;
    ASSERT_GE( N.at(0).at(0), 0.999999997-eps );
    ASSERT_LE( N.at(0).at(0), 0.999999997+eps );
    ASSERT_GE( N.at(1).at(0), 0.000069344-eps );
    ASSERT_LE( N.at(1).at(0), 0.000069344+eps );
    ASSERT_GE( N.at(2).at(0), 0.000030063-eps );
    ASSERT_LE( N.at(2).at(0), 0.000030063+eps );
    ASSERT_GE( N.at(0).at(1),-0.000069344-eps );
    ASSERT_LE( N.at(0).at(1),-0.000069344+eps );
    ASSERT_GE( N.at(1).at(1), 0.999999975-eps );
    ASSERT_LE( N.at(1).at(1), 0.999999975+eps );
    ASSERT_GE( N.at(2).at(1), 0.000011906-eps );
    ASSERT_LE( N.at(2).at(1), 0.000011906+eps );
    ASSERT_GE( N.at(0).at(2),-0.000030062-eps );
    ASSERT_LE( N.at(0).at(2),-0.000030062+eps );
    ASSERT_GE( N.at(1).at(2),-0.000011908-eps );
    ASSERT_LE( N.at(1).at(2),-0.000011908+eps );
    ASSERT_GE( N.at(2).at(2), 0.999999999-eps );
    ASSERT_LE( N.at(2).at(2), 0.999999999+eps );
}

TEST( TestMatrices, NutationTest2 )
{
    Time t( 2005, 4, 19, 0, 0, 0 );
    MatrixD N = nutation( t );
    double eps = 0.0001;
    ASSERT_GE( N.at(0).at(0), 0.999999999-eps );
    ASSERT_LE( N.at(0).at(0), 0.999999999+eps );
    ASSERT_GE( N.at(1).at(0),-0.000031216-eps );
    ASSERT_LE( N.at(1).at(0),-0.000031216+eps );
    ASSERT_GE( N.at(2).at(0),-0.000013928-eps );
    ASSERT_LE( N.at(2).at(0),-0.000013928+eps );
    ASSERT_GE( N.at(0).at(1),-0.000032122-eps );
    ASSERT_LE( N.at(0).at(1),-0.000032122+eps );
    ASSERT_GE( N.at(1).at(1), 0.999999999-eps );
    ASSERT_LE( N.at(1).at(1), 0.999999999+eps );
    ASSERT_GE( N.at(2).at(1), 0.000042516-eps );
    ASSERT_LE( N.at(2).at(1), 0.000042516+eps );
    ASSERT_GE( N.at(0).at(2), 0.000013926-eps );
    ASSERT_LE( N.at(0).at(2), 0.000013926+eps );
    ASSERT_GE( N.at(1).at(2),-0.000042516-eps );
    ASSERT_LE( N.at(1).at(2),-0.000042516+eps );
    ASSERT_GE( N.at(2).at(2), 0.999999999-eps );
    ASSERT_LE( N.at(2).at(2), 0.999999999+eps );
}

TEST( TestMatrices, NutationTest3 )
{
    Time t( 1992, 4, 19, 0, 0, 0 );
    MatrixD N = nutation( t );
    double eps = 0.0001;
    ASSERT_GE( N.at(0).at(0), 0.999999999-eps );
    ASSERT_LE( N.at(0).at(0), 0.999999999+eps );
    ASSERT_GE( N.at(1).at(0), 0.000071226-eps );
    ASSERT_LE( N.at(1).at(0), 0.000071226+eps );
    ASSERT_GE( N.at(2).at(0), 0.000030882-eps );
    ASSERT_LE( N.at(2).at(0), 0.000030882+eps );
    ASSERT_GE( N.at(0).at(1),-0.000071226-eps );
    ASSERT_LE( N.at(0).at(1),-0.000071226+eps );
    ASSERT_GE( N.at(1).at(1), 0.999999999-eps );
    ASSERT_LE( N.at(1).at(1), 0.999999999+eps );
    ASSERT_GE( N.at(2).at(1), 0.000048261-eps );
    ASSERT_LE( N.at(2).at(1), 0.000048261+eps );
    ASSERT_GE( N.at(0).at(2),-0.000030882-eps );
    ASSERT_LE( N.at(0).at(2),-0.000030882+eps );
    ASSERT_GE( N.at(1).at(2),-0.000004828-eps );
    ASSERT_LE( N.at(1).at(2),-0.000004828+eps );
    ASSERT_GE( N.at(2).at(2), 0.999999999-eps );
    ASSERT_LE( N.at(2).at(2), 0.999999999+eps );
}

TEST( TestMatrices, PolarTest1 )
{
    Time t( 2005, 4, 19, 0, 0, 0 );
    MatrixD W = polar( t );
    double eps = 0.005;
    Vector v1 = { 0.99999999999997, 0.000000000011734, 0.0000002540569133 };
    Vector v2 = {-0.00000000001207, 0.999999999999131, 0.0000013181841582 };
    Vector v3 = {-0.00000025405691,-0.000001318184158, 0.999999999999099  };
    MatrixD expW = { v1, v2, v3 };

    for( int i=0; i<3; ++i )
        for( int j=0; j<3; ++j )
        {
            ASSERT_GE( W.at(i).at(j), expW.at(i).at(j)-eps );
            ASSERT_LE( W.at(i).at(j), expW.at(i).at(j)+eps );
        }
}

TEST( TestMatrices, RotationTest1 )
{
    Time t( 2005, 4, 19, 0, 0, 0 );
    MatrixD G = rotation( t );
    double eps = 0.005;

    Vector v1 = { -0.889486850514503, 0.456960767201945, 0 };
    Vector v2 = { -0.456960767201945, -0.889486850514503, 0 };
    Vector v3 = { 0, 0, 1 };
    MatrixD expG = { v1, v2, v3 };

    for( int i=0; i<3; ++i )
        for( int j=0; j<3; ++j )
        {
            EXPECT_GE( G.at(i).at(j), expG.at(i).at(j)-eps );
            EXPECT_LE( G.at(i).at(j), expG.at(i).at(j)+eps );
        }
}

TEST( TestMatrices, RotationTest2 )
{
    Time t( 2010, 5, 5, 0, 0, 0 );
    MatrixD G = rotation( t );
    double eps = 0.005;

    Vector v1 = { -0.734183129027552, 0.678951495359803, 0 };
    Vector v2 = { -0.678951495359803, -0.734183129027552, 0 };
    Vector v3 = { 0, 0, 1 };
    MatrixD expG = { v1, v2, v3 };

    for( int i=0; i<3; ++i )
        for( int j=0; j<3; ++j )
        {
            EXPECT_GE( G.at(i).at(j), expG.at(i).at(j)-eps );
            EXPECT_LE( G.at(i).at(j), expG.at(i).at(j)+eps );
        }
}

// http://hpiers.obspm.fr/eop-pc/index.php?index=matrice&lang=en
TEST( TestMatrices, ConversionMatrixTest1 )
{
    Time t( 2005, 4, 19, 0, 0, 0 );
    MatrixD M = conversionMatrix( t );
    double eps = 0.005;
    Vector v1 = {-0.890012675378, 0.455935507093, 0.000501038234 };
    Vector v2 = {-0.455935545703,-0.890012795690, 0.000040898181 };
    Vector v3 = { 0.000464577372,-0.000192041241, 0.999999873644 };
    MatrixD expM = { v1, v2, v3 };
    for( int i=0; i<3; ++i )
        for( int j=0; j<3; ++j )
        {
            EXPECT_GE( M.at(i).at(j), expM.at(i).at(j)-eps );
            EXPECT_LE( M.at(i).at(j), expM.at(i).at(j)+eps );
        }
}

TEST( TestMatrices, ConversionMatrixTest2 )
{
    Time t( 2010, 5, 5, 0, 0, 0 );
    MatrixD M = conversionMatrix( t );
    double eps = 0.005;
    Vector v1 = {-0.735797525627, 0.677200803753, 0.001035701345 };
    Vector v2 = {-0.677201160082,-0.735797926542, 0.000008992792 };
    Vector v3 = { 0.000768156828,-0.000694761278, 0.999999463621 };
    MatrixD expM = { v1, v2, v3 };
    for( int i=0; i<3; ++i )
        for( int j=0; j<3; ++j )
        {
            EXPECT_GE( M.at(i).at(j), expM.at(i).at(j)-eps );
            EXPECT_LE( M.at(i).at(j), expM.at(i).at(j)+eps );
        }
}

int main( int argc, char **argv )
{
    testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();
}
