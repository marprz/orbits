#ifndef GSINTEGRATOR_H
#define GSINTEGRATOR_H

#include <list> // moze niepotrzebne
#include <deque>
#include <iomanip>
//#include "Integration.h"
//#include "Acceleration.h"
#include "Def.h"
#include "Conversion.cpp"
#include "Taylor.cpp"
#include "Integrator.h"

typedef std::vector< std::vector< Vector > > InitialVectors;

class GJIntegrator : public Integrator {
  public:
//    GJIntegrator();
    GJIntegrator( Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)( int, double, std::vector< Vector >) );
    GJIntegrator( const Vector& r0, const Vector& v0, Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)(int, double, std::vector< Vector >) );
//    Vector Acceleration( const Vector& position );
    InitialVectors InitializeVectorsGJ( int nb, double h );
    Vector Calcs0( const Vector& v0, const double& h );
    Vector CalcS0( const Vector& r0, const double& h );
    Vector Calcsn( const int& index );
    Vector CalcNexts( const Vector& prev_s, const Vector& prev_a, const Vector& curr_a );
    Vector CalcSn( const int& index );
    Vector CalcNextS( const Vector& prev_S, const Vector& prev_s, const Vector& prev_a );
    std::vector< Vector > CalcSumb();
    std::vector< Vector > CalcSuma();
    Vector Calcb5();
    Vector Calca5();
    void printData( std::string name1 = "positions", std::string name2 = "velocities", std::string name3 = "accelerations" );
    bool IsAccelerationConverged(); // TODO temporary !!!!!!!!
    void Startup( const double& h, std::vector< Vector >& corrected_rs, std::vector< Vector >& corrected_vs, std::vector< Vector >& updated_as );
    Vector PredictR( const double& h );
    Vector PredictV( const double& h );
    Vector MidCorrectR( const int& i, const double& h );
    Vector CorrectR( const double& h );
    Vector CorrectV( const double& h );
    void InitializeGaussJacksonCoefficients();
    void Algorithm( int points, double h );

    Vector (*Acceleration)(Vector);
    std::vector< Vector > (*Taylor)( int, double, std::vector< Vector > );
    std::vector< Vector > (*TaylorV)( int, double, std::vector< Vector > );

    MatrixD a, b;
    std::deque< Vector > sn, Sn; 
    double ra; // apogee [km]
    double rp; // perigee [km]
    double sa; // semi-major axis
    double T; // period ~12h
    bool isInitKnown;
    int current_n = 4;
};

void GJIntegrator::InitializeGaussJacksonCoefficients()
{
    b = { { 19087.0/89600, -427487.0/725760, 3498217.0/3628800, -500327.0/403200, 6467.0/5670, -2616161.0/3628800, 24019.0/80640, -263077.0/3628800, 8183.0/1036800 },
        {8183.0/1036800, 57251.0/403200, -1106377.0/3628800, 218483.0/725760, -69.0/280, 530177.0/3628800, -210359.0/3628800, 5533.0/403200, -425.0/290304},
        {-425.0/290304, 76453.0/3628800, 5143.0/57600, -660127.0/3628800, 661.0/5670, -4997.0/80640, 83927.0/3628800, -19109.0/3628800, 7.0/12800},
        {7.0/12800, -23173.0/3628800, 29579.0/725760, 2497.0/57600, -2563.0/22680, 172993.0/3628800, -6463.0/403200, 2497.0/725760, -2497.0/7257600},
        {-2497.0/7257600, 1469.0/403200, -68119.0/3628800, 252769.0/3628800, 0, -252769.0/3628800, 68119.0/3628800, -1469.0/403200, 2497.0/7257600},
        {2497.0/7257600, -2497.0/725760, 6463.0/403200, -172993.0/3628800, 2563.0/22680, -2497.0/57600, -29579.0/725760, 23173.0/3628800, -7.0/12800},
        {-7.0/12800, 19109.0/3628800, -83927.0/3628800, 4997.0/80640, -661.0/5670, 660127.0/3628800, -5143.0/57600, -76453.0/3628800, 425.0/290304},
        {425.0/290304, -5533.0/403200, 210359.0/3628800, -530177.0/3628800, 69.0/280, -218483.0/725760, 1106377.0/3628800, -57251.0/403200, -8183.0/1036800},
        {-8183.0/1036800, 263077.0/3628800, -24019.0/80640, 261616.0/3628800, -6467.0/5670, 500327.0/403200, -3498217.0/3628800, 427487.0/725760, -19087.0/89600},
        {25713.0/89600, -9401029.0/3628800, 5393233.0/518400, -9839609.0/403200, 167287.0/4536, -135352319.0/3628800, 10219841.0/403200, -40987771.0/3628800, 3288521.0/1036800 }};

    a = { { 3250433.0/53222400, 572741.0/5702400, -8701681.0/39916800, 4026311.0/13305600, -917039.0/3193344, 7370669.0/39916800, -1025779.0/13305600, 754331.0/39916800, -330157.0/159667200},
        {-330157.0/159667200, 530113.0/6652800, 518887.0/19958400, -27631.0/623700, 44773.0/1064448, -531521.0/19958400, 109343.0/9979200, -1261.0/475200, 45911.0/159667200},
        {45911.0/159667200, -185839.0/39916800, 171137.0/1900800, 73643.0/39916800, -25775.0/3193344, 77597.0/13305600, -98911.0/39916800, 24173.0/39916800, -3499.0/53222400},
        {-3499.0/53222400, 4387.0/4989600, -35039.0/4989600, 90817.0/950400, -20561.0/3193344, 2117.0/9979200, 2059.0/6652800, -317.0/2851200, 317.0/22809600},
        {317.0/22809600, -2539.0/13305600, 55067.0/39916800, -326911.0/39916800, 14797.0/152064, -326911.0/39916800, 55067.0/39916800, -2539.0/13305600, 317.0/22809600},
        {317.0/22809600, -317.0/2851200, 2059.0/6652800, 2117.0/9979200, -20561.0/3193344, 90817.0/950400, -35039.0/4989600, 4387.0/4989600, -3499.0/53222400},
        {-3499.0/53222400, 24173.0/39916800, -98911.0/39916800, 77597.0/13305600, -25775.0/3193344, 73643.0/39916800, 171137.0/1900800, -185839.0/39916800, 45911.0/159667200},
        {45911.0/159667200, -1261.0/475200, 109343.0/9979200, -531521.0/19958400, 44773.0/1064448, -27631.0/623700, 518887.0/19958400, 530113.0/6652800, -330157.0/159667200},
        {-330157.0/159667200, 754331.0/39916800, -1025779.0/13305600, 7370669.0/39916800, -917039.0/3193344, 4026311.0/13305600, -8701681.0/39916800, 572741.0/5702400, 3250433.0/53222400},
        {3250433.0/53222400, -11011481.0/19958400, 6322573.0/2851200, -8660609.0/1663200, 25162927.0/3193344, -159314453.0/19958400, 18071351.0/3326400, -24115843.0/9979200, 103798439.0/159667200 }};
}
#endif // GSINTEGRATOR_H
