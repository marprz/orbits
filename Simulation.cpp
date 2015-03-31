#include <list> // moze niepotrzebne
#include <deque>
#include <iomanip>
#include "Integration.h"
#include "Acceleration.h"
#include "Conversion.cpp"
#include "Taylor.cpp"

typedef std::vector< std::vector< Vector > > InitialVectors;

InitialVectors InitializeVectorsGJ( int nb, double h )
{
    InitialVectors init;
    double ra = kRE+20394; // apogee [km]
    double rp = kRE+19970; // perigee [km]
    double sa = (rp+ra)/2;  // semi-major axis
    double T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };
    Vector a0 = Acceleration( r0 );

    std::vector< Vector > fr = { r0, v0, a0 };
    std::vector< Vector > fv = { v0, a0 };
    // TWO-BODY TAYLOR:
    std::vector< Vector > rs = taylorSeriesTwoBody( 10, h, fr ); // potrzeba 9 pozycji do algorytmu. 10-ta jest do obliczenia predkosci jako pochodnych pozycji
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
Vector Calcs0( const Vector& v0, const MatrixD& b, const std::vector< Vector >& as, const double& h )
{
    Vector ret;
    Vector sum = { 0, 0, 0 };
    for( int i=0; i<as.size(); ++i )
        sum = sum + b.at(4).at(i)*as.at(i);
    ret = v0/h-sum;
    print( ret, "s0" );
    return ret;
}

// equation 86
Vector CalcS0( const Vector& r0, const MatrixD& a, const std::vector< Vector >&as, const double& h )
{
    Vector ret;
    Vector sum = { 0, 0, 0 };
    for( int i=0; i<as.size(); ++i )
        sum = sum + a.at(4).at(i)*as.at(i);
    ret = r0/pow(h,2)-sum;
    return ret; 
}
// equation 75
void Calcsn( std::deque< Vector >& s, const std::vector< Vector >& a, const int& index )
{
    if( index == 0 )
        std::cout << "ERROR! Calcsn - index should not be 0! " << std::endl;
    if( index < 0 )
    {
        // example: index=-1 => a(n) = 3, a(n+1) = 4; 3=index+4, 4=index+5
        Vector sn = s.front() - ( a.at(index+5) + a.at(index+4) )/2;
        s.push_front( sn );
    } 
    else
    {
        Vector sn = s.back() + ( a.at(index+4-1) + a.at(index+5-1) )/2;
        s.push_back( sn );
    }
}

// equation 86
void CalcSn( const std::deque< Vector >& s, std::deque< Vector >& S, const std::vector< Vector >& a, const int& index )
{
    std::cout << "function CalcSn: index = " << index  << std::endl;
    if( index == 0 )
        std::cout << "ERROR! CalcSn - index should not be 0! " << std::endl;
    if( index < 0 )
    {
        // example: index=-1 => a(n) = 3, a(n+1) = 4; 3=index+4, 4=index+5
        Vector Sn = S.front() - s.at(index+5) +(a.at(index+5))/2;
        S.push_front( Sn );
    } 
    else
    {
        std::cout << "s.size() = " << s.size() << ", a.size() = " << a.size() << std::endl;
        Vector Sn = S.back() + s.at(index+4) + a.at(index+4)/2;
        std::cout << "new Sn calculated" << std::endl;
        S.push_back( Sn );
        print( S.back(), "new calculated last Sn ");
    }
}

// table 5
std::vector< Vector > CalcSumb( const MatrixD& b, std::vector< Vector > as )
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
std::vector< Vector > CalcSuma( const MatrixD& a, std::vector< Vector > as )
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
Vector Calcb5( const MatrixD& b, std::vector< Vector > as )
{
    Vector ret = { 0, 0, 0 };
    int as_size = as.size(); 
    for( int i=0; i<b.at(9).size(); ++i )
        ret = ret + b.at(9).at(9-i-1)*as.at(as_size-i-1);
    return ret;
}

// PREDICT - step 6.
Vector Calca5( const MatrixD& a, std::vector< Vector > as )
{
    Vector ret = { 0, 0, 0}; 
    int as_size = as.size(); 
    for( int i=0; i<a.at(9).size(); ++i )
        ret = ret + a.at(9).at(9-i-1)*as.at(as_size-i-1) ;
    return ret;
}

void SimulationFromZeroGaussJackson()
{    int n = 9;
    double h = 60; // step size in seconds
    h = 5*h; // step is 5 min
    InitialVectors init = InitializeVectorsGJ( n, h );
    std::vector< Vector > rs = init.at(0);
    std::vector< Vector > vs = init.at(1);
    std::vector< Vector > as = init.at(2);

    print( rs, "two-body positions");
    /*print( vs, "two-body velocities");
    print( as, "two-body accelerations");*/
    std::vector< Vector >::iterator it;
    std::cout << "positions: ";
    for( it = rs.begin(); it != rs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << "velocities: ";
    for( it = vs.begin(); it != vs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << "accelerations: ";
    for( it = as.begin(); it != as.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;

    MatrixD a, b;
    GaussJacksonCoefficients( a, b );

    std::deque< Vector > sn, Sn; // TODO powinna byc lepsza struktura, zamiast tworzenia nowego wektora z listy

    bool is_acceleration_converged = false;
    int counter = 0; 
    int pointNb = 200; 
    while( pointNb > 0 )
    {
        --pointNb;
        while( !is_acceleration_converged )
        {
            std::vector< Vector > sum_bn, sum_an;

            // sn:
            Vector s0 = Calcs0( vs.at(4), b, as, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
            sn.push_back( s0 );
            for( int i=1; i<5; ++i )
            {
                Calcsn( sn, as, -i ); 
                Calcsn( sn, as, i ); 
            }

            // Sn:
            Vector S0 = CalcS0( rs.at(4), a, as, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
            Sn.push_back( S0 );
            for( int i=1; i<5; ++i )
            {
                CalcSn( sn, Sn, as, -i ); // CalcSn( const std::vector< Vector >& s, std::list< Vector >& S, const std::vector< Vector >& a, const int& index )
                CalcSn( sn, Sn, as, i ); 
            }

            // step 3b  iii.
            sum_bn = CalcSumb( b, as ); 
            sum_an = CalcSuma( a, as );

            // step 3b iv 
            std::vector< Vector > corrected_vs; // eq. 74
            std::vector< Vector > corrected_rs; // eq. 87
            std::vector< Vector > updated_as; 
            for( int i=0; i<9; ++i )
            {   
                if( i == 4 )
                {
                    corrected_vs.push_back( vs.at(4) );
                    corrected_rs.push_back( rs.at(4) );
                }
                else 
                {
                    corrected_vs.push_back( h*(sn.at(i)+sum_bn.at(i)) );
                    corrected_rs.push_back( h*h*(Sn.at(i)+sum_an.at(i)) );
                }
                updated_as.push_back( Acceleration( corrected_rs.back() ) );
            }

            std::cout << std::setprecision( 5 );
    /*        print( rs, "two-body positions");
            print( corrected_rs, "corrected_rs" );
            print( vs, "two-body velocities");
            print( corrected_vs, "corrected_vs" );
            print( as, "two-body accelerations");
            print( updated_as, "updated_as" );
            std::cout << "positions: ";
            for( it = corrected_rs.begin(); it != corrected_rs.end(); ++it )
                std::cout << norm(*it) << " ";
            std::cout << std::endl;
            std::cout << "velocities: ";
            for( it = corrected_vs.begin(); it != corrected_vs.end(); ++it )
                std::cout << norm(*it) << " ";
            std::cout << std::endl;
            std::cout << "updated_accelerations: ";
            for( it = updated_as.begin(); it != updated_as.end(); ++it )
                std::cout << norm(*it) << " ";
            std::cout << std::endl;*/

            /* najlepsze jest rozwiniecie Taylor'a. Dlaczego? TODO 
            rs = corrected_rs;
            vs = corrected_vs;*/
            ++counter;
            if( counter > 0 )
                is_acceleration_converged = true; // TODO temporary
        }

        //***************************************************************** PREDICT 
        //
        for( int c = 0; c<1; ++c )
        //for( int c = 0; c<2; ++c )
        {
            // 4. Calculate S(n+1)
            Vector S = Sn.back() + sn.back() + as.back()/2;
            Sn.push_back( S );
            Sn.pop_front();

            // 5. & 6. Calculate b5 and a5
            Vector b5 = Calca5( b, as );
            Vector a5 = Calca5( a, as );

            // 7. Calculate pred_vn, pred_rn
            Vector pred_vn = { 0, 0, 0 }; // eq. 77
            Vector pred_rn = { 0, 0, 0 }; // eq. 88
            pred_vn = h*(sn.back() + as.back() + b5 );
            pred_rn = h*h*(S + a5 );

            print( pred_vn, "predicted v");
            print( pred_rn, "predicted r");
    
            vs.push_back( pred_vn );
            rs.push_back( pred_rn );

            // EVALUATE - CORRECT
            // 8. Evaluate the acceleration
            as.push_back( Acceleration( pred_rn ) );
    
            int index = 5;
            bool v_and_r_converged = false;
            int it_number = 1;
            int max_it_number = 5;
            // 10. while v and r have not converged 
            //     and the maximum nb of corrector iterations is not exceeded:
            while( !v_and_r_converged && it_number <= max_it_number )
            {
                std::cout << "in while " << std::endl;
                // TODO Calcsn zapisuje nowe sn na koniec kolejki
                Calcsn( sn, as, index ) ; 
                sn.pop_front(); // TODO temporary prawdopodobnie 
                std::cout << "calcSn done" << std::endl;
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
                // step 10 c.
                // pominiete, bo w 10b policzone dla wszyskich k
                // step 10 d. 
                // sum: sn.back(), sum_b4, sum_a4, 10c to obtain:
                // v (79) and r (89)
                std::cout << "size as = " << as.size() << std::endl;
                Vector next_v = h*(sn.back() + sum_b4);
                Vector next_r = h*h*(Sn.back() + sum_a4);
        //        Vector new_rn = h*(sn.back() + )
            
                rs.push_back( next_r );
                vs.push_back( next_v );

                ++it_number;
                v_and_r_converged = true;
            } // while max_it_number and ...  
        }
    } // pointNb = 0
    std::cout << "positions: ";
    for( it = rs.begin(); it != rs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
    std::cout << "velocities: ";
    for( it = vs.begin(); it != vs.end(); ++it )
        std::cout << norm(*it) << " ";
    std::cout << std::endl;
        std::cout << "updated_accelerations: ";
        for( it = as.begin(); it != as.end(); ++it )
            std::cout << norm(*it) << " ";
        std::cout << std::endl;


        std::vector< std::vector< Vector > > positionsVector;
        positionsVector.push_back( rs );
        saveVectorToFile( positionsVector, "positionsVector" );
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
    saveVectorToFile( positions, "positions");
*/
}

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
}

int main()
{
    SimulationFromZeroGaussJackson();
//    SimulationFromZeroNystrom();
//    SimulationNGA();  // nie dziala - V0 nie jest znane
//    SimulationSin();
    return 0;
}
