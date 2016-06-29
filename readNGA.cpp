#include <iostream>
#include <fstream>
#include <string>
#include "Def.h"
#include "Conversion.cpp"
#include "Acceleration.h"

std::pair< Vector, Vector > readNGA( std::string in_file_name, std::string out_file_name )
{
    int n = 1;
    std::vector< std::vector< Vector > > transformed_positions;
    transformed_positions.resize( n );

    std::fstream file;
    file.open( in_file_name.c_str(), std::ios::in );

    std::cout << "file " << in_file_name.c_str() << std::endl;
    Vector r0 = { 0, 0, 0 };
    Vector v0 = { 0, 0, 0 };
    if( file.is_open() )
    {
        int year, m, d, h, mt, id;
        std::string type, temp;
        double sec, x, y, z;
        double temp2;
        while( file >> temp >> year >> m >> d >> h >> mt >> sec )
        {
            std::cout << year << " " << m << " " << d << " " << h << " " << mt << " " << sec << std::endl;
            Time t( year, m, d, h, mt, sec );
           
            for( int j=0; j<64; ++j )
            {
                file >> type >> id >> x >> y >> z >> temp2;
                if( type == std::string("V") && id<=n )
                //if( type == std::string("P") )
                {
                    Vector pos = convert( t, x, y, z );
                    (transformed_positions.at(id-1)).push_back( pos );
                }
                else
                    if( id == 1 && h==mt==sec==0 && type == std::string("V") )
                    {
                        v0 = convert( t,x/1000,y/1000,z/1000 );
                    }
            }
        }
    
        std::cout << "transformed_positions(0).size() = " << transformed_positions.at(0).size() << std::endl;
        file.close();
        //saveVectorToFile( transformed_positions, out_file_name );
        for( int i=0; i<n; ++i )
        {
            out_file_name = "sat";
            out_file_name += std::to_string(i+1);
            saveVectorToFile( transformed_positions.at(i), out_file_name );
        }
    }
    else 
        std::cout << "ERROR! File " << in_file_name << " IS NOT open" << std::endl;

    r0 = transformed_positions.at(0).at(0);
    return std::make_pair( r0, v0 );
}

int main( int argc, char* argv[] )
{
    std::string out_file_name;
    if( argc < 2 )
        out_file_name = std::string("out_file");
    else
        out_file_name = argv[1];
    std::string fileName;
    if( argc < 3 ) 
        fileName = std::string( "NGA17702" );
    else
        fileName = argv[2];
    std::pair< Vector, Vector > initRV = readNGA( fileName, out_file_name );
//    print( initRV.first, "r0" );
//    print( initRV.second, "v0" );
    return 0;
}
