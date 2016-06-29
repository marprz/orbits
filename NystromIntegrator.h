#ifndef NYSTROM_INTEGRATOR_H
#define NYSTROM_INTEGRATOR_H
#include "Integrator.h"

class NystromIntegrator : public Integrator
{
  public:
    NystromIntegrator();
    void Startup( const Vector& position, const Vector& velocity );
    Vector Acceleration( const Vector& position );
    void InitializeParameters();
    void Algorithm( int points, double h );
    
    // positions, velocities and accelerations:
    std::vector< Vector > rs;
    std::vector< Vector > vs;
    std::vector< Vector > as;

    // parameters:
    Vector ALPHA;
    std::vector< Vector > A;
    Vector AL; 
    double ra; // apogee [km]
    double rp; // perigee [km]
    double sa;  // semi-major axis
    double T; // period ~12h
};
#endif// NYSTROM_INTEGRATOR_H
