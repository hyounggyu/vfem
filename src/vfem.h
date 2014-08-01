#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

const double ns           = 1.45; // ns=silica index, must be always maximum
const double a_bc_marker  = 1.0; //  Artificial Bcs Marker -fixed
const double s2_bc_marker = 2.0; //  Symmetry Bcs-2 Marker -fixed
const double s3_bc_marker = 3.0; //  Symmetry Bcs-2 Marker -fixed
const double de_bc_marker = 4.0; //  Di-eletric Bcs Marker -fixed

const double PI            = 3.14159265;
const int    MAX_FNAME_LEN = 200;

enum POLAR { X_POLAR, Y_POLAR };
