#include "GJIntegrator.cpp"

int main()
{
    std::cout << "TEST INTEGRATOR" << std::endl;
    GJIntegrator gs_integrator;
    double h = 10;
    int points = 1000;
    gs_integrator.Algorithm( points, h );
    return 0;
}

