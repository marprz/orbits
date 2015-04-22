#include "AdamsIntegrator.h"

AdamsIntegrator::AdamsIntegrator( Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)( int, double, std::vector< Vector >) )
{
    ra = kRE+20394; // apogee [km]
    rp = kRE+19970; // perigee [km]
    sa = (rp+ra)/2;  // semi-major axis
    T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    Taylor = p_Taylor;
    TaylorV = p_TaylorV;
    Acceleration = p_Acceleration;
    isInitKnown = false;
}

AdamsIntegrator::AdamsIntegrator( const Vector& r0, const Vector& v0, Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)(int, double, std::vector< Vector >) )
{
    ra = kRE+20394; // apogee [km]
    rp = kRE+19970; // perigee [km]
    sa = (rp+ra)/2;  // semi-major axis
    T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h
    setPos0( r0 );
    setVel0( v0 );
    isInitKnown = true;
    
    Taylor = p_Taylor;
    TaylorV = p_TaylorV;
    Acceleration = p_Acceleration;
}


InitialVectors AdamsIntegrator::InitializeVectorsGJ( int nb, double h )
{
    InitialVectors init;

    // initial position, velocity and acceleration:
    if( !isInitKnown )
    {
        double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
        setPos0( { rp, 0, 0 } );
        setVel0( { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 } );
    }

    print( r0, "r0" );
    print( Acceleration( r0 ) , "acceleration(r0)" );
    std::cout << "fr: " << std::endl;
    std::vector< Vector > fr = { r0, v0, Acceleration( r0 ) };
    std::vector< Vector > fv = { v0, Acceleration( r0 ) };

    std::cout << "Taylor: " << std::endl;
    // TWO-BODY TAYLOR:
    rs = Taylor( 10, h, fr ); // potrzeba 9 pozycji do algorytmu. 10-ta jest do obliczenia predkosci jako pochodnych pozycji
    vs = TaylorV( 9, h, fr );

    std::vector< Vector >::iterator it;
    /*for( int i=0; i<rs.size()-1; ++i )
    {
        vs.push_back( 1/h*(rs.at(i+1)-rs.at(i)) );
    }*/
    rs.pop_back(); // ostatni element byl potrzebny tylko do liczenia predkosci
    for( it = rs.begin(); it != rs.end(); ++it )
    {
        as.push_back( Acceleration( *it ) );
    }
    init.push_back( rs );
    init.push_back( vs );
    init.push_back( as );
    return init;

}

void AdamsIntegrator::Startup( const double& h, std::vector< Vector >& corrected_rs, std::vector< Vector >& corrected_vs, std::vector< Vector >& updated_as )
{

}

void AdamsIntegrator::Algorithm( int points, double h )
{    
    int n = 9;
    double h1; // step size in seconds
    h1 = h;
    double h2 = h1;
    InitialVectors init = InitializeVectorsGJ( n, h1 );
    rs = init.at(0);
    vs = init.at(1);
    as = init.at(2);

    for( int i=0; i<points; ++i )
    {
        int curr_ind = vs.size();
        Vector next_v = vs.at(curr_ind-1) + h*( as.at(curr_ind-1)*251/720 + as.at(curr_ind-2)*646/720 - as.at(curr_ind-3)*264/720 + as.at(curr_ind-4)*106/720 - as.at(curr_ind-5)*19/720);
        vs.push_back( next_v );
        Vector next_r = rs.at(curr_ind-1) + h*( vs.at(curr_ind-1)*251/720 + vs.at(curr_ind-2)*646/720 - vs.at(curr_ind-3)*264/720 + vs.at(curr_ind-4)*106/720 - vs.at(curr_ind-5)*19/720);
        rs.push_back( next_r );
        as.push_back( Acceleration( next_r ) );
    }

    saveVectorToFile( rs, "positions" );
}


