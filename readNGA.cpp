#include <iostream>
#include <fstream>
#include <string>
#include "Def.h"
#include "Conversion.cpp"

void saveVectorToFile( const std::vector< std::vector< Vector3 > >& positions, const std::string& outName )
{
    int satellites_nb = positions.size();
    int positions_nb = positions.at(0).size();
    std::cout << "satellites_nb = " << satellites_nb << ", positions_nb = " << positions_nb << std::endl;
    std::fstream outFile;
    outFile.open( outName.c_str(), std::ios::out | std::ios::in );
    outFile << "function p = " << outName << "\n";
    outFile << "p=zeros(" << satellites_nb << ", " << positions_nb << ", 3 );\n";
    for( int i=0; i<satellites_nb; ++i )
    {
        int temp_pos_nb = positions.at(i).size();
        if( temp_pos_nb > 0 )
        {
            std::cout << "i=" << i << ", temp_pos_nb=" << temp_pos_nb << std::endl;
            outFile << "p(" << i+1 << ",:,:) = [ ";
            for( int j=0; j<temp_pos_nb-1; ++j )
            {
                double x = ((positions.at(i)).at(j)).at(0);
                double y = ((positions.at(i)).at(j)).at(1);
                double z = ((positions.at(i)).at(j)).at(2);
                outFile << x << " " << y << " " << z << ";\n";
            }
            outFile << ((positions.at(i)).at(temp_pos_nb-1)).at(0) << " " << ((positions.at(i)).at(temp_pos_nb-1)).at(1) << " " << ((positions.at(i)).at(temp_pos_nb-1)).at(2) << "];\n";
        }
    }
    outFile.close();
}

int main( int argc, char* argv[] )
{
    std::string out_file_name;
    if( argc < 2 )
        out_file_name = std::string("out_file");
    else
        out_file_name = argv[1];

    std::cout << "file: " << out_file_name << std::endl;
    std::vector< std::vector< Vector3 > > transformed_positions;
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
                    Vector3 pos = convert( t, x, y, z );
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
