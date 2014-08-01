#include "vfem.h"
#include "matrix.h"

Matrix::Matrix () {

}

int Matrix::resize ( int size ) {
    elements.resize(size);
    return 0;
}

double Matrix::find ( int row, int col ) {
    int val = 0;
    std::list<ELEMENT>::iterator it = elements[col].begin();

    while ( it != elements[col].end() && it->row < row )
        ++it;

    if (it->row == row)
        return it->val;
    else
        return val;
}

int Matrix::set_element ( int row, int col, double val ) {
    ELEMENT tmp;
    std::list<ELEMENT>::iterator it = elements[col].begin();
    
    while ( it != elements[col].end() && it->row < row ) {
        ++it;
    }
    
    if ( it->row == row ) {
        it->val = val;
    } else {
        tmp.row = row, tmp.val = val;
        elements[col].insert(it, tmp);
    }
    
    return 0;
}

int Matrix::get_rowsize ( int col ) {
    return elements[col].size();
}

int Matrix::get_row_val ( int col, int *irow, double *val ) {
    std::list<ELEMENT>::iterator it = elements[col].begin();

    int i;
    for ( i=0; it != elements[col].end(); ++i, ++it ) {
        irow[i] = it->row;
        val[i] = it->val;
    }

    return 0;
}

#ifdef DEBUG
int Matrix::print ( ) {
    std::list<ELEMENT>::iterator it;
    
    for ( int i=0; i < elements.size(); i++) {
        it = elements[i].begin();
        while( it != elements[i].end() ) {
            printf("%d\t%d\t%f\n", it->row, i, it->val);
            ++it;
        }
    }
}
#endif
