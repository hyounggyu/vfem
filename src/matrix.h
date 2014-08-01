#include <vector>
#include <list>

class Matrix {
    typedef struct {
        int row;
        double val;
    } ELEMENT;
    std::vector< std::list<ELEMENT> > elements;

public:
    Matrix ();
    int resize( int );
    double find ( int, int );
    int set_element ( int, int, double );
    int get_rowsize ( int );
    int get_row_val ( int, int *, double *);
#ifdef DEBUG
    int print ();
#endif
};
