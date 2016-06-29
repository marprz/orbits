#include "NystromIntegrator.h"

NystromIntegrator::NystromIntegrator()
{
    InitializeParameters();
    ra = kRE+20394; // apogee [km]
    rp = kRE+19970; // perigee [km]
    sa = (rp+ra)/2;  // semi-major axis
    T = 2*M_PI*sqrt( pow(sa,3)/kmuE ); // period ~12h

    double Vp = sqrt( kmuE*(2/rp-1/sa)); // speed at perigee
    std::cout << "perigee: " << rp << std::endl;
    Vector r0 = { rp, 0, 0 };
    Vector v0 = { 0, Vp*sqrt(2)/2, Vp*sqrt(2)/2 };

    Startup( r0, v0 );
}

void NystromIntegrator::Startup( const Vector& position, const Vector& velocity )
{
    rs.push_back( position );
    vs.push_back( velocity );
    as.push_back( Acceleration( position ) );
}

Vector NystromIntegrator::Acceleration( const Vector& position )
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
    return acc;
}

void NystromIntegrator::Algorithm( int points, double h )
{
    for( int i=0; i<points; ++i )
    {
        Vector r0 = rs.back();
        Vector v0 = vs.back();
        Vector tempV = ALPHA.at(0)*h*v0;
        tempV = r0+tempV;
//    Vector temp = { -1*r0 };
        Vector temp = Acceleration( tempV );

        const int s = 12; // from Hairer
        std::array< Vector, s > k;
        k.at(0) = temp;
        for( int i=1; i<s-1; ++i )
        {
            Vector sum = { 0, 0, 0 };
            for( int j=0; j<i-1; ++j )
            {
                sum = sum + (A.at(i-1).at(j))*(k.at(j));
            }
            Vector par = pow(h,2)*sum + r0+ALPHA.at(i)*h*v0;
            k.at(i) = Acceleration( par );
    //        k.at(i) = -1*par;
        }

        Vector sum2 = { 0, 0, 0 };
        for( int i=0; i<s-1; ++i )
        {
            sum2 = sum2 + A.at(s-1).at(i)*k.at(i);
        }
        Vector r1 = h*h*sum2 + r0+h*v0;

        Vector sum3 = { 0, 0, 0 };
        for( int i=0; i<s-1; ++i )
        {
            sum3 = sum3 + AL.at(i)*k.at(i);
        }
        Vector v1 = h*sum3+v0;

        rs.push_back( r1 );
        vs.push_back( v1 );
        as.push_back( Acceleration( r1 ) );
    }

    std::cout << "obliczono " << rs.size() << " pozycji " << std::endl;
    std::vector< std::vector< Vector > > positions;
    positions.push_back( rs );
    saveVectorToFile( positions, "positionsNystrom");

}

void NystromIntegrator::InitializeParameters()
{
    ALPHA.resize(11);
    ALPHA.at(0) = 0;
    ALPHA.at(1) = 8.48093651182655889791*0.01;
    ALPHA.at(2) = 1.69618730236531177958*0.1;
    ALPHA.at(3) = 8.82527661964732346426*0.1;
    ALPHA.at(4) = 3.57384241759677451843*0.1;
    ALPHA.at(5) = 1.17472338035267653574*0.1;
    ALPHA.at(6) = 6.42615758240322548157*0.1;
    ALPHA.at(7) = ALPHA.at(3);
    ALPHA.at(8) = ALPHA.at(2);
    ALPHA.at(9) = 0;
    ALPHA.at(10) = 1;

    A.resize(12);
    for( int i=0; i<12; ++i )
    {
        A.at(i).resize(11);
        for( int j=0; j<11; ++j )
            A.at(i).at(j) = 0;
    }
    
    A.at(1).at(0)  = 3.59631420588164200966*0.001;
    A.at(2).at(0)  = 4.79508560784218934622*0.001;
    A.at(2).at(1)  = 9.59017121568437869244*0.001;
    A.at(3).at(0)  = 1.87733546050627049772;
    A.at(3).at(1)  =-4.32661231970111510045;
    A.at(3).at(2)  = 2.83870439626131304487;
    A.at(4).at(0)  = 1.94709467473823596644*0.01;
    A.at(4).at(1)  = 0;
    A.at(4).at(2)  = 4.42810873502123510231*0.01;
    A.at(4).at(3)  = 1.09714031475079632769*0.0001;
    A.at(5).at(0)  = 4.64206697543819328333*0.001;
    A.at(5).at(1)  = 0;
    A.at(5).at(2)  = 2.93579155497804175774*0.001;
    A.at(5).at(3)  = 2.76435808114904080285*0.00001;
    A.at(5).at(4)  =-7.05627009491629727298*0.0001;
    A.at(6).at(0)  =-6.96029489417689828096*0.001;
    A.at(6).at(1)  = 0;
    A.at(6).at(2)  =-1.90054512658393594723*0.1;
    A.at(6).at(3)  = 7.68088148461498565180*0.1;
    A.at(6).at(4)  = 1.18704223034434343286*0.1;
    A.at(6).at(5)  = 2.84020002739066989630*0.1;
    A.at(7).at(0)  = 8.26380420350248571364*0.01;
    A.at(7).at(1)  = 0;
    A.at(7).at(2)  = 3.60887996724131292687*0.1;
    A.at(7).at(3)  = 1.01208078941331907108*0.001;
    A.at(7).at(4)  = 6.79462118331253280025*0.01;
    A.at(7).at(5)  =-2.04227957300564960765*0.1;
    A.at(7).at(6)  = 8.11711629853386060148*0.01;
    A.at(8).at(0)  = 2.38468700370378669757*0.01;
    A.at(8).at(1)  = 9.59017121568437869244*0.001;
    A.at(8).at(2)  = 1.52655867995367678583*0.1;
    A.at(8).at(3)  =-5.41908820084510828926*0.001;
    A.at(8).at(4)  =-3.93767256187270659404*0.01;
    A.at(8).at(5)  =-1.40338104799298795297*0.1; 
    A.at(8).at(6)  = 1.00094410252962054139*0.01;
    A.at(8).at(7)  = 3.41682516901140790025*0.001;
    A.at(9).at(0) =-4.25990467506340701948*0.01;
    A.at(9).at(1) = 0;
    A.at(9).at(2) =-5.93485416695720137147*0.1;
    A.at(9).at(3) = 9.81481906263695126914*0.001;
    A.at(9).at(4) = 7.13421080958931560696*0.01;
    A.at(9).at(5) = 2.91368066505515770072*0.1;
    A.at(9).at(6) =-1.57502473657776550315*0.01;
    A.at(9).at(7) =-6.97794897116397988241*0.001;
    A.at(9).at(8) = 2.86287666119249964845*0.1;
    A.at(10).at(0) =-6.14180937931931543943*0.1;
    A.at(10).at(1) = 0;
    A.at(10).at(2) =-1.75100055986614758156;
    A.at(10).at(3) =-6.76603738760857112582*0.001;
    A.at(10).at(4) = 3.95581030796451903257*0.1;
    A.at(10).at(5) = 1.13603914070927980965;
    A.at(10).at(6) = 5.64165276602374526099*0.01;
    A.at(10).at(7) = 3.58989223779823259146*0.01;
    A.at(10).at(8) = 7.48011913641736205396*0.1;
    A.at(10).at(9)= 5*0.1;
    A.at(11).at(0) = 0;
    A.at(11).at(1) = 0;
    A.at(11).at(2) =-8.30381269763468822042*0.01;
    A.at(11).at(3) = 0;
    A.at(11).at(4) = 1.78280368337326917339*0.1;
    A.at(11).at(5) = 1.67007309146871573764*0.1;
    A.at(11).at(6) = 9.91488201804162591694*0.01;
    A.at(11).at(7) = 2.22301690020519163947*0.01;
    A.at(11).at(8) = 8.30381269763468822042*0.01;
    A.at(11).at(9)= 3.33333333333333333333*0.01;
    A.at(11).at(10)= 0;
    
    AL.resize(11);
    AL.at(0) = 0;
    AL.at(1) = 0;
    AL.at(2) =-1*0.1;
    AL.at(3) = 0;
    AL.at(4) = 2.77429188517743176508*0.1;
    AL.at(5) = 1.89237478148923490158*0.1;
    AL.at(6) = 2.77429188517743176508*0.1;
    AL.at(7) = 1.89237478148923490158*0.1;
    AL.at(8) = 1*0.1;
    AL.at(9)= 3.33333333333333333333*0.01;
    AL.at(10)= 3.33333333333333333333*0.01;
}


