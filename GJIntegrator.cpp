#include "GJIntegrator.h"

GJIntegrator::GJIntegrator()
{
    ra = kRE+20394; // apogee [km]
    rp = kRE+19970; // perigee [km]
    sa = (rp+ra)/2;  // semi-major axis
    T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h
}

Vector GJIntegrator::Acceleration( const Vector& position )
{
    Vector acc = { 0, 0, 0 };
    double x = position.at(0);
    double y = position.at(1);
    double z = position.at(2);
    double r_mag = sqrt( x*x + y*y + z*z );
    double nuR3 = -kmuE/(pow(r_mag,3));
    double pp = 3/2 *kmuE*kJ2*kRE2/(pow(r_mag,5));
    double z2r2 = 5*pow(position.at(2)/r_mag,2);
    Vector c = { 1, 1, 3 };
    for( int i=0; i<3; ++i )
    {
        acc.at(i) = nuR3*position.at(i)-pp*position.at(i)*(c.at(i)-z2r2);
    }
    print( acc, "ACCELERATION");
    return acc;
}

InitialVectors GJIntegrator::InitializeVectorsGJ( int nb, double h )
{
    InitialVectors init;
    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };
    Vector a0 = Acceleration( r0 );

    std::vector< Vector > fr = { r0, v0, a0 };
    std::vector< Vector > fv = { v0, a0 };
    // TWO-BODY TAYLOR:
    std::vector< Vector > rs = taylorSeriesTwoBody( 10, h, fr ); // potrzeba 9 pozycji do algorytmu. 10-ta jest do obliczenia predkosci jako pochodnych pozycji
//    std::vector< Vector > rs = taylorSeriesSin( 10, h ); // potrzeba 9 pozycji do algorytmu. 10-ta jest do obliczenia predkosci jako pochodnych pozycji
    std::vector< Vector > vs;
    std::vector< Vector > as;
    std::vector< Vector >::iterator it;
    for( int i=0; i<rs.size()-1; ++i )
    {
        vs.push_back( 1/h*(rs.at(i+1)-rs.at(i)) );
    }
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

// equation 75 s0 = C1'
Vector GJIntegrator::Calcs0( const Vector& v0, const double& h )
{
    bool debug = true;
    Vector ret;
    Vector sum = { 0, 0, 0 };
    std::cout << "Calculating s0:" << std::endl;
    for( int i=0; i<as.size(); ++i )
    {
        sum = sum + b.at(4).at(i)*as.at(i);
        if( debug )
        {
            std::cout << "b.at(4).at(" << i << ") = " << b.at(4).at(i) << std::endl << "as.at(" << i << ")=";
            print( as.at(i), "");
            print( sum, "sum");
        }
    }
    ret = v0/h-sum;
    print( ret, "s0" );
    return ret;
}

// equation 86
Vector GJIntegrator::CalcS0( const Vector& r0, const double& h )
{
    bool debug = true;
    Vector ret;
    Vector sum = { 0, 0, 0 };
    std::cout << "Calculating S0: " << std::endl;
    for( int i=0; i<as.size(); ++i )
    {
        sum = sum + a.at(4).at(i)*as.at(i);        
        if( debug )
        {
            std::cout << "a.at(4).at(" << i << ") = " << a.at(4).at(i) << std::endl << "as.at("<< i << ")=";
            print( as.at(i), "");
            print( sum, "sum");
        }
    }

    ret = r0/pow(h,2)-sum;
    return ret; 
}

// equation 75
Vector GJIntegrator::Calcsn( const int& index )
{
    Vector s;
    if( index == 0 )
        std::cout << "ERROR! Calcsn - index should not be 0! " << std::endl;
    if( index < 0 )
    {
        // example: index=-1 => a(n) = 3, a(n+1) = 4; 3=index+4, 4=index+5
        s = sn.front() - ( as.at(index+5) + as.at(index+4) )/2;
    } 
    else
    {
        s = sn.back() + ( as.at(index+4-1) + as.at(index+5-1) )/2;
    }
    return s;
}

Vector GJIntegrator::CalcNexts( const Vector& prev_s, const Vector& prev_a, const Vector& curr_a )
{
    Vector s = prev_s + (prev_a + curr_a)/2;
    return s;
}

// equation 86 
Vector GJIntegrator::CalcSn( const int& index )
{
    Vector S;
    if( index == 0 )
        std::cout << "ERROR! CalcSn - index should not be 0! " << std::endl;
    if( index < 0 )
    {
        // example: index=-1 => a(n) = 3, a(n+1) = 4; 3=index+4, 4=index+5
        S = Sn.front() - sn.at(index+5) +(as.at(index+5))/2;
    } 
    else
    {
        S = Sn.back() + sn.at(index+4) + as.at(index+4)/2;
    }
    return S;
}

Vector GJIntegrator::CalcNextS( const Vector& prev_S, const Vector& prev_s, const Vector& prev_a )
{
    Vector S = prev_S + prev_s + prev_a/2;
    return S;
}

// table 5
std::vector< Vector > GJIntegrator::CalcSumb()
{
    std::vector< Vector > ret;
    for( int i=0; i<10; ++i ) // n=-4,...,4
    {
        Vector sum = { 0, 0, 0 };
        for( int k=0; k<9; ++k )
        {
            sum = sum + b.at(i).at(k)*as.at(k);
        }
        ret.push_back( sum );
    }
    return ret;
}

// table 6
std::vector< Vector > GJIntegrator::CalcSuma()
{
    std::vector< Vector > ret;
    for( int i=0; i<9; ++i ) // n=-4,...,4
    {
        Vector sum = { 0, 0, 0 };
        for( int k=0; k<9; ++k )
        {
            sum = sum + a.at(i).at(k)*as.at(k);
        }
        ret.push_back( sum );
    }
/*    std::cout << "CalcSuma: " << std::endl;
    for( int i=0; i<ret.size(); ++i )
        print( ret.at(i), "ret.at(i)" );
    std::cout << std::endl;*/
    return ret;
}

// PREDICT - step 5.
Vector GJIntegrator::Calcb5()
{
    Vector ret = { 0, 0, 0 };
    int as_size = as.size(); 
    for( int i=0; i<b.at(9).size(); ++i )
        ret = ret + b.at(9).at(9-i-1)*as.at(as_size-i-1);
    return ret;
}

// PREDICT - step 6.
Vector GJIntegrator::Calca5()
{
    Vector ret = { 0, 0, 0}; 
    int as_size = as.size(); 
    for( int i=0; i<a.at(9).size(); ++i )
    {
        ret = ret + a.at(9).at(9-i-1)*as.at(as_size-i-1) ;
    }
    return ret;
}

void GJIntegrator::printData( std::string name1, std::string name2, std::string name3 )
{
    print( rs, name1);
    print( vs, name2);
    print( as, name3);
    std::vector< Vector >::const_iterator it;
    std::cout << name1 << ": ";
    for( it = rs.begin(); it != rs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << name2 << ": ";
    for( it = vs.begin(); it != vs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << name3 << ": ";
    for( it = as.begin(); it != as.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
}

bool GJIntegrator::IsAccelerationConverged() // TODO temporary !!!!!!!!
{
    return true; // FAKE TODO zrobic tu prawdziwa funkcje!
}

void GJIntegrator::Startup( const double& h, std::vector< Vector >& corrected_rs, std::vector< Vector >& corrected_vs, std::vector< Vector >& updated_as, Vector r0, Vector v0 )
{
    // sn:
    Vector s0 = Calcs0( v0, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
    sn.push_back( s0 );
    for( int i=1; i<5; ++i )
    {
         sn.push_front( Calcsn( -i ) ); 
         sn.push_back( Calcsn( i ) ); 
    }

    // Sn:
    Vector S0 = CalcS0( r0, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
    Sn.push_back( S0 );
    for( int i=1; i<5; ++i )
    {
        Sn.push_front( CalcSn( -i ) ); // CalcSn( const std::vector< Vector >& s, std::list< Vector >& S, const std::vector< Vector >& a, const int& index )
        Sn.push_back( CalcSn( i ) ); 
    }

    // step 3b  iii.
    std::vector< Vector > sum_bn, sum_an;
    sum_bn = CalcSumb(); 
    sum_an = CalcSuma();

    updated_as.clear();
    // step 3b iv 
    for( int i=0; i<9; ++i )
    {   
        if( i == 4 )
        {
            corrected_vs.push_back( v0 );
            corrected_rs.push_back( r0 );
        }
        else 
        {
            corrected_vs.push_back( h*(sn.at(i)+sum_bn.at(i)) );
            corrected_rs.push_back( h*h*(Sn.at(i)+sum_an.at(i)) );
        }
        updated_as.push_back( Acceleration( corrected_rs.back() ) );
    }
}

Vector GJIntegrator::PredictR( const double& h )
{
    // 6. Calculate a5
    Vector a5 = Calca5();

    // 7. Calculate pred_vn, pred_rn
    Vector pred_rn = h*h*(Sn.back() + a5 );

//    print( pred_vn, "predicted v");
    print( pred_rn, "predicted r");
    return pred_rn;
}

Vector GJIntegrator::PredictV( const double& h )
{
    // 5. Calculate b5 
    Vector b5 = Calca5();

    // 7. Calculate pred_vn, pred_rn
    Vector pred_vn = h*(sn.back() + as.back() + b5 );

    print( pred_vn, "predicted v");
    return pred_vn;
}

Vector GJIntegrator::CorrectR( const double& h )
{
    Vector corr_rn;
    std::cout << " EVALUATE - CORRECT" << std::endl;
    int index = as.size()-5;
    Vector sum_a4 = { 0, 0, 0 };
    for( int k=0; k<9; ++k )
    //for( int k=0; k<8; ++k )
    {
        // step 10 b.
        sum_a4 = sum_a4 + a.at(8).at(k)*as.at(index-4+k);
    }
    corr_rn = h*h*(Sn.back() + sum_a4);
    print( corr_rn, "corrected r");
           
    as.pop_back();
    as.push_back( Acceleration( corr_rn ) );

    return corr_rn;
}

Vector GJIntegrator::CorrectV( const double& h )
{
    Vector corr_vn;
    std::cout << " EVALUATE - CORRECT" << std::endl;
    int index = 5;
    int s_as = as.size();
    Vector s = CalcNexts( sn.back(), as.at( s_as-2 ), as.back()  ) ; 
    sn.push_back( s ) ; 
    Vector sum_b4 = { 0, 0, 0 };
    Vector sum_a4 = { 0, 0, 0 };
//                if( it_number == 1 )
    {
        for( int k=0; k<9; ++k )
        //for( int k=0; k<8; ++k )
        {
            // step 10 b.
            sum_b4 = sum_b4 + b.at(8).at(k)*as.at(index-4+k);
            sum_a4 = sum_a4 + a.at(8).at(k)*as.at(index-4+k);
        }
    }
    corr_vn = h*(sn.back() + sum_b4);
    print( corr_vn, "corrected v");
            
//                vs.push_back( next_v );

    return corr_vn;
}

void GJIntegrator::Algorithm( int points, double h )
{
//    points = static_cast< int >( T/h );
    int n = 9;
    double h1 = 5; // step size in seconds
    h1 = h;
    double h2 = h1;
    int pointNb1 = points; 
    int pointNb2 = pointNb1;;
    InitialVectors init = InitializeVectorsGJ( n, h1 );
    rs = init.at(0);
    vs = init.at(1);
    as = init.at(2);

    InitializeGaussJacksonCoefficients();

    std::cout << "po inicjalizacji " << std::endl;

    bool is_acceleration_converged = false;
    int counter = 1; 
    int tempCounter = 0; 
    
    while( !is_acceleration_converged )
    {
        std::vector< Vector > sum_bn, sum_an;
        std::cout << std::setprecision( 5 );
        std::vector< Vector > corrected_vs; // eq. 74
        std::vector< Vector > corrected_rs; // eq. 87
        std::vector< Vector > updated_as; 

        Startup( h1, corrected_rs, corrected_vs, updated_as, rs.at(4), vs.at(4) );
        if( tempCounter == 0 ) 
        {
            is_acceleration_converged = IsAccelerationConverged(); // TODO temporary !!!!!!!!
        }
        ++tempCounter;
        for( int i=0; i<9; ++i )
        {
            Vector tempDiff = as.at(i) - updated_as.at(i);
            print( tempDiff, "difference(i)" );
        }
        for( int i=0; i<9; ++i )
        {
            Vector tempDiff = rs.at(i) - corrected_rs.at(i);
            print( tempDiff, "differenceR" );
        }

        rs = corrected_rs;
        vs = corrected_vs; 
        for( int i = 0; i < rs.size(); ++i )
        {
            as.at(i) = Acceleration( rs.at(i) );
        }
    }

    MatrixD::iterator iterV;
    for( iterV = rs.begin(); iterV != rs.end(); ++iterV )
        print( *iterV, "position(i)");
    for( iterV = as.begin(); iterV != as.end(); ++iterV )
        print( *iterV, "acceleration(i)");



    int iterPoints = points;
    while( iterPoints > 0 )
    {
        --iterPoints;

        //*********************************************************** PREDICT 
        //
        // 4. Calculate S(n+1)
        Vector S = Sn.back() + sn.back() + as.back()/2;
        Sn.push_back( S );
        Sn.pop_front();

        Vector pred_vn = { 0, 0, 0 }; // eq. 77
        Vector pred_rn = { 0, 0, 0 }; // eq. 88
        // c - liczba iteracji
        pred_rn = PredictR( h1 );
        pred_vn = PredictV( h1 );

//        Predict( Sn, sn, a, b, as, h, pred_rn, pred_vn );
        rs.push_back( pred_rn );
        vs.push_back( pred_vn );

        // EVALUATED ACCELERATION: 
        print( Acceleration( pred_rn ), " last_as" );
        as.push_back( Acceleration( pred_rn ) );

        for( int k = 0; k<4; ++k )
        {

            Vector corr_vn; // eq. 79
            Vector corr_rn; // eq. 89
            int max_it_number = 4;
            corr_rn = CorrectR( h1 ); 
            //corr_rn = Correct( counter, as, pred_rn, corr_vn, max_it_number, Sn, sn, h, a, b ); 
            Vector pos_diff = corr_rn-pred_rn;
            print( pos_diff, "difference in positions");

       //     vs.push_back( corr_vn );
            rs.pop_back();
            rs.push_back( corr_rn );
        }

        int sn_size = sn.size();
        Vector next_sn = sn.back() + ( as.back() + as.at(sn_size-2) )/2;
        sn.push_back( next_sn );
        sn.pop_front();
    } // pointNb = 0
//    printData( rs, vs, as );

    std::vector< std::vector< Vector > > positionsVector;
    std::vector< std::vector< Vector > > accVector;
    positionsVector.push_back( rs );
    accVector.push_back( as );
    std::cout << "rs.size = " << rs.size();
    std::cout << "positionsVector.size() = " << positionsVector.size() << std::endl;
    std::cout << "positionsVector.at(0).size() = " << positionsVector.at(0).size() << std::endl;
    std::cout << "positionsVector.at(0).at(0).size() = " << positionsVector.at(0).at(0).size() << std::endl;
    saveVectorToFile( positionsVector, "positionsVector" );
    saveVectorToFile( accVector, "accVector" );
/*    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positions");*/
}

    /*
void SimulationFromZeroNystrom()
{
    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };

    // vectors of positions, velocities
    std::vector< Vector > rs;
    std::vector< Vector > vs;

    rs.push_back( r0 );
    vs.push_back( v0 );

    double h = 60; // step size in seconds
    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positionsNystrom");

}

// r - position
// v - velocity 
// a - acceleration
void SimulationNGA() // nie dziala
{
    // vectors of positions, velocities, accelerations
    std::vector< Vector > rs;
    std::vector< Vector > vs;
    std::vector< Vector > as;

    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h
//    double Vp = sqrt( muE*(2/rp-1/sa)); // speed at perigee

    Time t( 2013, 12, 10, 0, 0, 0 );
    Vector r0_TRS = { -15210.576086, 2042.480101, 21649.168938 };
    Vector r0 = convert( t, r0_TRS );
    Vector v0_TRS = { -7773.910720, -26449.272652, -2851.405296 };
    Vector v0 = 0.001*v0_TRS ;
    //Vector v0 = convert( t, v0_TRS );

    double r = norm(r0);
    std::cout << "Vp = " << sqrt( kmuE*(2/r-1/sa) ) << std::endl;
    std::cout << "v0 = " << v0.at(0) << " " << v0.at(1) << " " << v0.at(2) << std::endl;
    rs.push_back( r0 );
    vs.push_back( v0 );
    as.push_back( Acceleration( r0 ) );

    double h = 60; // step size: 5 minutes
    for( int ti = h; ti<T; ti += h )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positionsNystrom");
}

void SimulationSin()
{
    std::vector< Vector > rs;
    std::vector< Vector > vs;
  
    Vector r0 = { 0 };
    Vector v0 = { 1 };
    rs.push_back( r0 );
    vs.push_back( v0 );
    double h = 0.02;
    for( int i=0; i<300; ++i )
    {
        std::vector< Vector > solved_ode = OdeNystrom( r0, v0, h );
        r0 = solved_ode.at(0);
        v0 = solved_ode.at(1);
       // a0 = solved_ode.at(2);
        rs.push_back( r0 );
        vs.push_back( v0 );
        //as.push_back( Acceleration( r0 ) );
    }
    std::cout << "wektor pozycji dla sinusa (h=0.1): ";
    for( int i=0; i<rs.size(); ++i )
    {
        std::cout << rs.at(i).at(0) << " (" << sin(i*h) << "), ";
    }
    std::cout << std::endl;
}*/
/*
int main()
{
    SimulationFromZeroGaussJackson();
//    SimulationFromZeroNystrom();
//    SimulationNGA();  // nie dziala - V0 nie jest znane
//    SimulationSin();
    return 0;
}*/
