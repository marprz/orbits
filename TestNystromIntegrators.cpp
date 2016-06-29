#include "NystromIntegrator.cpp"

void NystromTest()
{
    NystromIntegrator nystrom_integrator;
    double h=10;
    int points = 10000;
    nystrom_integrator.Algorithm( points, h );
}

int main()
{
    NystromTest();
    std::cout << "TEST NYSTROM INTEGRATOR" << std::endl;
    return 0;
}

