#include <iostream>
#include <fstream>
#include <string>
#include "Def.h"
#include "Conversion.cpp"
#include "Acceleration.h"

int main( int argc, char* argv[] )
{
    std::string out_file_name;
    if( argc < 2 )
        out_file_name = std::string("out_file");
    else
        out_file_name = argv[1];

    std::cout << "file: " << out_file_name << std::endl;
    std::vector< std::vector< Vector > > transformed_positions;
    transformed_positions.resize( 32 );

    std::fstream file;
    std::string fileName( "NGA17702" );
    file.open( fileName.c_str(), std::ios::in );
    if( file.is_open() )
    {
        int year, m, d, h, mt, id;
        std::string type;
        double sec, x, y, z, temp;
        while( file >> year >> m >> d >> h >> mt >> sec )
        {
            Time t( year, m, d, h, mt, sec );
//            t.print();
            for( int j=0; j<62; ++j )
            {
    //            std::cout << "j=" << j << std::endl;
                file >> type >> id >> x >> y >> z;
    //            std::cout << type << " " << id << " " << x << " " << y << " " << z << std::endl;
                if( type == std::string("P") )
                {
                    Vector pos = convert( t, x, y, z );
                    (transformed_positions.at(id-1)).push_back( pos );
                }
            }
        }
    
        saveVectorToFile( transformed_positions, out_file_name );
        file.close();
    }
    else 
        std::cout << "ERROR! File " << fileName << " IS NOT open" << std::endl;
    return 0;
}
