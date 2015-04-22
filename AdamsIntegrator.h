#ifndef ADAMSINTEGRATOR_H
#define ADAMSINTEGRATOR_H

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

class AdamsIntegrator : public Integrator {
  public:
//    GJIntegrator();
    AdamsIntegrator( Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)( int, double, std::vector< Vector >) );
    AdamsIntegrator( const Vector& r0, const Vector& v0, Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)(int, double, std::vector< Vector >) );
    InitialVectors InitializeVectorsGJ( int nb, double h );
    void printData( std::string name1 = "positions", std::string name2 = "velocities", std::string name3 = "accelerations" );
    void Startup( const double& h, std::vector< Vector >& corrected_rs, std::vector< Vector >& corrected_vs, std::vector< Vector >& updated_as );
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
#endif // ADAMSINTEGRATOR_H
