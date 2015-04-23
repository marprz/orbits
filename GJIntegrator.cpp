#include "GJIntegrator.h"

GJIntegrator::GJIntegrator( Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)( int, double, std::vector< Vector >) )
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

GJIntegrator::GJIntegrator( const Vector& r0, const Vector& v0, Vector (*p_Acceleration)(Vector), std::vector< Vector > (*p_Taylor)( int, double, std::vector< Vector > ), std::vector< Vector > (*p_TaylorV)(int, double, std::vector< Vector >) )
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

InitialVectors GJIntegrator::InitializeVectorsGJ( int nb, double h )
{
    std::cout << "h = " << h << std::endl;
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

// equation 75 s0 = C1'
Vector GJIntegrator::Calcs0( const Vector& v0, const double& h )
{
    std::cout << "Calcs0\n as.size() = " << as.size() << std::endl;
    print( v0, "v0" );
    bool debug = true;
    Vector ret;
    Vector sum = { 0, 0, 0 };
    for( int i=0; i<as.size(); ++i )
    {
        sum = sum + b.at(4).at(i)*as.at(i);
    }
    ret = v0/h-sum;
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
        std::cout << "Calculating sn(" << index << ") - takes as.at(" << index+5 << ") and as.at(" << index+4 << ")" << std::endl;
        s = sn.front() - ( as.at(index+5) + as.at(index+4) )/2;
    } 
    else
    {
        std::cout << "Calculating sn(" << index << ") - takes as.at(" << index+4-1 << ") and as.at(" << index+5-1 << ")" << std::endl;
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
        std::cout << "Calculating Sn(" << index << ") - takes sn.at(" << index+5 << ") and as.at(" << index+5 << ")" << std::endl;
        S = Sn.front() - sn.at(index+5) +(as.at(index+5))/2;
    } 
    else
    {
        std::cout << "Calculating Sn(" << index << ") - takes sn.at(" << index+4 << ") and as.at(" << index+4 << ")" << std::endl;
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

void GJIntegrator::Startup( const double& h, std::vector< Vector >& corrected_rs, std::vector< Vector >& corrected_vs, std::vector< Vector >& updated_as )
{
    std::cout << "h = " << h << std::endl;
    std::cout << "Begin Startup " << std::endl;
    // sn:
    sn.clear();
    Vector s0 = Calcs0( v0, h ); //Calcs0( Vector v0, MatrixD b, std::vector< Vector > as, double h )
    sn.push_back( s0 );
    for( int i=1; i<5; ++i )
    {
         sn.push_front( Calcsn( -i ) ); 
         sn.push_back( Calcsn( i ) ); 
    }

    // Sn:
    Sn.clear();
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
            Vector suma = { 0, 0, 0 };
            for( int j=0; j<9; ++j )
            {
                suma = suma + (a.at(i)).at(j)*as.at(j);
            }
            corrected_vs.push_back( h*(sn.at(i)+sum_bn.at(i)) ); // eq. 74 - mid-corrector
            corrected_rs.push_back( h*h*(Sn.at(i)+ suma ) ); // eq. 87 - mid-corrector
            //corrected_rs.push_back( h*h*(Sn.at(i)+sum_an.at(i)) ); // eq. 87 - mid-corrector
        }
        updated_as.push_back( Acceleration( corrected_rs.back() ) );
    }
    std::cout << "End of Startup" << std::endl << std::endl;
}

Vector GJIntegrator::PredictR( const double& h )
{
    std::cout << "h = " << h << std::endl;
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

Vector GJIntegrator::MidCorrectR( const int& i, const double& h )
{
    Vector ak = a.at(i); 
    Vector _sum = { 0, 0, 0 };
    for( int j=0; j<9; ++j )
        _sum = _sum + ak.at(j)*as.at(i);
    return ( h*h*(Sn.at(i)) + _sum );
}

Vector GJIntegrator::CorrectR( const double& h )
{
    Vector corr_rn;
    int as_size = as.size();
    Vector sum_a4 = { 0, 0, 0 };
    for( int k=0; k<9; ++k )
    //for( int k=0; k<8; ++k )
    {
        // step 10 b.
        sum_a4 = sum_a4 + (a.at(8)).at(k)*as.at(as_size-9+k);
    }
    corr_rn = h*h*(Sn.back() + sum_a4);
    
    return corr_rn;
}

Vector GJIntegrator::CorrectV( const double& h )
{
    Vector corr_vn;
    int index = 5;
    int s_as = as.size();
    Vector s = CalcNexts( sn.back(), as.at( s_as-2 ), as.back()  ) ; 
    sn.push_back( s ) ; 
    Vector sum_b4 = { 0, 0, 0 };
//                if( it_number == 1 )
    {
        for( int k=0; k<9; ++k )
        //for( int k=0; k<8; ++k )
        {
            // step 10 b.
            sum_b4 = sum_b4 + b.at(8).at(k)*as.at(index-4+k);
        }
    }
    corr_vn = h*(sn.back() + sum_b4);
    return corr_vn;
}

void GJIntegrator::Algorithm( int points, double h )
{
    points = T/h+1;
    std::cout << "points: " << points << std::endl;
//    points = static_cast< int >( T/h );
    int n = 9;
    int pointNb1 = points; 
    int pointNb2 = pointNb1;;
    std::cout << "InitializeVectorsGJ: " << std::endl;
    InitialVectors init = InitializeVectorsGJ( n, h );
    std::cout << "InitializeVectorsGJ - done " << std::endl;
    rs = init.at(0);
    vs = init.at(1);
    as = init.at(2);

    InitializeGaussJacksonCoefficients();

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

        Startup( h, corrected_rs, corrected_vs, updated_as );
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
        
        for( int i = 0; i < rs.size(); ++i )
        {
            as.at(i) = Acceleration( rs.at(i) );
        }


 
        std::cout << "rs.size() = " << rs.size() << std::endl << std::endl;
/*        std::cout << "Compare positions with mid-corrector: " << std::endl;
        for( int i=0; i<9; ++i )
        {
            print( rs.at(i), "rs:" );
            print( newrs.at(i), "new" );
            std::cout << std::endl;
        }*/
    }

    std::cout << "After startup rs.size()=" << rs.size() << std::endl;
    MatrixD::iterator iterV;
    for( iterV = rs.begin(); iterV != rs.end(); ++iterV )
        print( *iterV, "position(i)");
    for( iterV = as.begin(); iterV != as.end(); ++iterV )
        print( *iterV, "acceleration(i)");
    
    std::vector< Vector > uncorrected_rs = rs;
    std::vector< Vector > uncorrected_vs = vs;
    int iterPoints = points;
    while( iterPoints > 0 )
    {
        --iterPoints;
        int pointInd = points-iterPoints; 
        std::cout << "\nCALCULATING POINT " << points-iterPoints << std::endl;

        //*********************************************************** PREDICT 
        //
        // 4. Calculate S(n+1)
        Vector S = Sn.back() + sn.back() + as.back()/2;
        Sn.push_back( S );
        Sn.pop_front();

        Vector pred_vn = { 0, 0, 0 }; // eq. 77
        Vector pred_rn = { 0, 0, 0 }; // eq. 88
        // c - liczba iteracji
        pred_rn = PredictR( h );
        pred_vn = PredictV( h );

        uncorrected_rs.push_back( pred_rn );
        uncorrected_vs.push_back( pred_vn );
//        Predict( Sn, sn, a, b, as, h, pred_rn, pred_vn );
        rs.push_back( pred_rn );
        vs.push_back( pred_vn );

        // EVALUATED ACCELERATION: 
        as.push_back( Acceleration( pred_rn ) );

        print( as.back(), "predicted as" );
        print( pred_rn, "predicted r" );
        
        double temph = (pointInd+4)*h; 
        Vector ex_r = { sin(temph), 0, 0 };
    //    print( ex_r, "exact r" );
        Vector corr_vn; // eq. 79
        Vector corr_rn; // eq. 89       

/*        // usunac przy odkomentowywaniu kodu ponizej
            int as_size = as.size();
            Vector next_sn = sn.back() + ( as.back() + as.at(as_size-2) )/2;
            sn.push_back( next_sn );
            print( sn.back(), "calculated sn" );
*/
        // Correction:
        for( int k = 0; k<1; ++k )
        {
            std::cout << "CORRECTION NB. " << k+1 << std::endl;

            if( k != 0 )
            {
                sn.pop_back();
            }
            int as_size = as.size();
            Vector next_sn = sn.back() + ( as.back() + as.at(as_size-2) )/2;
            sn.push_back( next_sn );
            print( sn.back(), "calculated sn" );

            corr_rn = CorrectR( h ); 
            corr_vn = CorrectV( h );

            Vector d_r = corr_rn - rs.back(); 
            print( d_r, "roznica w polozeniu: ");

            Vector d_a = as.back() - Acceleration( corr_rn );
            print( d_a, "przyspieszenie - roznica");
            print( Acceleration( corr_rn ), "calc acceleration" );

            as.pop_back();
            as.push_back( Acceleration( corr_rn ) );

            rs.pop_back();
            rs.push_back( corr_rn );
            vs.pop_back();
            vs.push_back( corr_vn );           
            
        }
//        print( ex_r, "exact r" );
        sn.pop_front();
    } // pointNb = 0
//    printData( rs, vs, as );

    std::cout << "h=" << h << std::endl;

    saveVectorToFile( rs, "rs", h );
//    saveVectorToFile( uncorrected_vs, "uncorrected_vs" );
//    saveVectorToFile( uncorrected_rs, "uncorrected_rs" );
    saveVectorToFile( vs, "vs" );
    saveVectorToFile( as, "as" );
}
