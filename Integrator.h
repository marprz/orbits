#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include "Def.h"

class Integrator
{
  public:
    void setPos0( const Vector& init_pos );
    void setVel0( const Vector& init_vel );
//    void printData();
  protected:
    Vector r0;
    Vector v0;
    Vector a0;
    MatrixD rs;
    MatrixD vs;
    MatrixD as;
};

void Integrator::setPos0( const Vector& init_pos )
{
    r0 = { init_pos };
}

void Integrator::setVel0( const Vector& init_vel )
{
    v0 = { init_vel };
}

#endif // INTEGRATOR_H
