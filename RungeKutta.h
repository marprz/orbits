#ifndef RK_INTEGRATOR_H
#define RK_INTEGRATOR_H
#include "Integrator.h"

class RungeKutta : public Integrator
{
  public:
    RungeKutta( Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)( int, double, std::vector< Vector >) );
    RungeKutta( const Vector& r0, const Vector& v0, Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)(int, double, std::vector< Vector >) );

    void Startup( const Vector& position, const Vector& velocity );
    Vector Acceleration( const Vector& position );
    void Algorithm( int points, double h );
    
    Vector (*Acceleration)(Vector);
    std::vector< Vector > (*Taylor)( int, double, std::vector< Vector > );
    std::vector< Vector > (*TaylorV)( int, double, std::vector< Vector > );


    // parameters:
    double ra; // apogee [km]
    double rp; // perigee [km]
    double sa;  // semi-major axis
    double T; // period ~12h
    bool isInitKnown;
};
#endif// RK_INTEGRATOR_H
