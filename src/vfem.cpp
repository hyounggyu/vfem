// ISOTROPIC/AnISOTROPIC VFEM with {Hx ,Hy, Hz}
// Ref 1: Performance comparsion of FE approaches for EM waveguides
//          by  S.Selleri  and M. Zaboli , 1997
// Vector Finite Element Method (Program);
// A Penalty function Method- BY S.SELLERI, 1997
// (For Real Effective Index problems)
//
// Khurram Naeem
// ONTL,GIST  
// S KOREA

// For calculating Caldding modes all Symmetry BCs combination in PCF.
// For Effective Mode area and Relative sensitivity of each guided mode

#include "yaml-cpp/yaml.h"
#include "vfem.h"

// k0 : wave-vector
int nmodes;
double wavelength, n1, n2, n3, k0, k02;
bool PEC_PEC;
POLAR polar;

int elements, **id, *elem_marker, nodes, *marker;
double **coord, *elem_area;

int read_input(const char *, std::string &, std::string &);
int read_mesh(std::string);
extern int run(void);
extern int post_process(std::string &);

int main (int argc,char *argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << *argv << " inputfile" << std::endl;
        exit(1);
    }

    std::string mesh_prefix, output;
    read_input(argv[1], mesh_prefix, output);
    read_mesh(mesh_prefix);

    run();

    post_process(output);

    return 0;
}

int read_input(const char *fname, std::string &mesh_prefix, std::string &output) {
    std::ifstream input;
    input.open(fname);

    if (!input.is_open()) {
        std::cerr << "Input file open error" << std::endl;
        exit(1);
    }

    YAML::Parser parser(input);
    YAML::Node doc;

    parser.GetNextDocument(doc);

    std::string polar_s, bc_s;

    doc["params"]["wavelength"]         >> wavelength;  // wavelength is in micro meter
    doc["params"]["n1"]                 >> n1;
    doc["params"]["clad_hole_index"]    >> n2;          // clad hole index -> 2nd maximum
    doc["params"]["core_hole_index"]    >> n3;          // core hole index
    doc["params"]["boundary_condition"] >> bc_s;
    doc["params"]["polarization"]       >> polar_s;
    doc["params"]["no_of_modes"]        >> nmodes;
    doc["files"]["mesh_prefix"]         >> mesh_prefix;
    doc["files"]["output"]              >> output;

    // set k0^2, k0 is wave vector
    k0  = 2*PI/wavelength;
    k02 = k0*k0;

    // <!> insert parameter check codes

    // set polarization
    if (polar_s.compare("X") == 0 || polar_s.compare("x") == 0) {
        polar = X_POLAR;
    } else if (polar_s.compare("Y") == 0 || polar_s.compare("y") == 0) {
        polar = Y_POLAR;
    } else {
        std::cerr << "Polarization input error: " << polar_s << std::endl;
        exit(1);
    }

    // set boundary condition
    // PEC-PMC combination will work if PEC_PEC==false.
    // else PEC-PEC or PMC-PMC Bcs applied.
    if (bc_s.compare("PEC-PEC") == 0 || bc_s.compare("PMC-PMC") == 0) {
        PEC_PEC = true;
    } else if (bc_s.compare("PEC-PMC") == 0) {
        PEC_PEC = false;
    } else {
        std::cerr << "Boundary condition input error: " << bc_s << std::endl;
    }

    return 0;
}

// Reading element file(Scanning Global node numbers)
int read_mesh(std::string mesh_prefix) {
    std::ifstream f;
    int tmp_int;
    float tmp_float;
    std::string fname, tmp_string;

    // 1. Reading element file(Scanning Global node numbers)
    fname = mesh_prefix + ".e";
    f.open(fname.c_str(), std::ifstream::in);
    if(!f.is_open()) {
        std::cerr << "Element file open error" << std::endl;
        exit(1);
    }
    
    // first line is number of elements
    f >> elements;

    elem_marker = new int[elements];
    elem_area   = new double[elements];
    id          = new int*[elements];
    for ( int i=0; i < elements; i++ )
        id[i] = new int[3];

    // <!> add line format comments
    for ( int i=0; i < elements; i++ ) {
        f >> tmp_string;
        f >> id[i][0] >> id[i][1] >> id[i][2];
        f >> tmp_int >> tmp_int >> tmp_int
          >> tmp_int >> tmp_int >> tmp_int;
        f >> tmp_float >> tmp_float;
        f >> tmp_int;
        elem_marker[i] = tmp_int;
    }

    f.close();

    // 2. Reading NODE file(Scanning coordinate of eash node) 
    fname = mesh_prefix + ".n";
    f.open(fname.c_str(), std::ifstream::in);
    if(!f.is_open()) {
        std::cerr << "Element file open error" << std::endl;
        exit(1);
    }

    // first line is number of nodes
    f >> nodes;

    marker = new int[nodes];

    coord = new double*[nodes];
    for ( int i=0; i < nodes; i++ )
        coord[i] = new double[2];
    
    for ( int i=0; i < nodes; i++ )
        f >> tmp_string >> coord[i][0] >> coord[i][1] >> marker[i];

    f.close();

    return 0;
}
