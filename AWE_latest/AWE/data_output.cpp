#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include "data_output.hpp"

using namespace std;


void data_output(int& iter, vector<SPH_particle>* particle_list)
{
    string file_num;
    string filename;
    file_num = to_string(iter);

    filename = "../../post_process/output_" + file_num + ".txt";

    fstream fs(filename, std::fstream::out);

    fs << "# Data output for iteration" << iter << "\n";
    fs << setw(20) << "Position_x";
    fs << setw(20) << "Position_y";
    fs << setw(20) << "Pressure";
    fs << setw(20) << "Velocity_x";
    fs << setw(20) << "Velocity_y";
    fs << setw(20) << "Density\n";

    for (auto p = particle_list->begin(); p != particle_list->end(); ++p)
    {
        fs << setw(20) << p->x[0];
        fs << setw(20) << p->x[1];
        fs << setw(20) << p->P;
        fs << setw(20) << p->v[0];
        fs << setw(20) << p->v[0];
        fs << setw(20) << p->rho << endl;
    }
    
    fs.flush();
    fs.close();
}
