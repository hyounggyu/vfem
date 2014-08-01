#include "areig.h"
#include "matrix.h"
#include "vfem.h"

int set_boundaries(void);
int set_matrix(int, int *, double *, double *);
int eigensystem(void);
int get_csc_matrix(Matrix &, int &, double *, int *, int *);

extern int nmodes;
extern double n1, n2, n3, k0, k02;
extern POLAR polar;
extern bool PEC_PEC;

extern int elements, **id, *elem_marker, nodes, *marker;
extern double **coord, *elem_area;

const int NNZ_PER_COL = 26;

int *bs1, *bs2, *bs3, *bc_marker, n, size;
double NxNx[3][3],NyNy[3][3],NxNy[3][3],NN[3][3],NNx[3][3],NNy[3][3];
Matrix S, T;

int nconv, nev;
double *EigValR, *EigValI, *EigVec;

int run() {
    // Element shape functions:
    // p[0] ,p[1] & p[2] are 1st, 2nd & 3rd global nodes no in Anticlockwise rotation of j-ELEMENT.
    int p[3];
    double x[3],y[3];

    S.resize(3*nodes);
    T.resize(3*nodes);

    set_boundaries();

    for ( int je=0; je < elements; je++ ) {

        for(int i=0; i < 3;i++) {
            p[i]=id[je][i];
            // coordinates of p[0], p[1] & p[2] of an element je are:
            x[i]=coord[p[i]][0];
            y[i]=coord[p[i]][1];
        }

        set_matrix(je, p, x, y);
    }

    eigensystem();

    return 0;
}

// Reading NODE file(Scanning coordinate of eash node)
// Noted that symetric Boundaries markers are fixed
// i.e bs2 = 2 for X-axis & bs3=3 for 2nd symm plane.
int set_boundaries(void) {
    bc_marker = new int[nodes];

    bs1 = new int[3*nodes];   // Artificial Bcs = 1
    bs2 = new int[3*nodes];   // Symmetry Bcs = 2
    bs3 = new int[3*nodes];   // Symmetry Bcs = 3

    // c1   : No. of unknown nodes except Artificial-1 BCs applied.
    // c12  : No. of unknown nodes except Artificial-1 BCs & Symmetry-2 BCs.
    // c13  : No. of unknown nodes except Artificial-1 BCs & Symmetry-3 BCs.
    // c123 : No. of unknown nodes except Artificial-1 BCs, Symmetry-2 BCs & Symmetry-3 BCs.
    int c1, c12, c13, c123, point0=1, ii=0;

    // un1  : For excluding  boundry (artificial 1) bcs
    // un2  : For excluding  boundry (symetry2 + artificial) bcs
    // un3  : For excluding  boundry (symetry3 + artificial) bcs
    // un23 : For excluding all boundry (symetry2 + symetry3 + artificial) bcs
    int un1=0, un2=0, un3=0, un23=0;

    for ( int i=0; i < nodes; i++ ) {
        if ((coord[i][0]==0)&&(coord[i][1]==0)&&(marker[i] == 3)) {
            if (PEC_PEC) {
                marker[i] = 1, point0 = i, ii = 1;
            } else if (polar == X_POLAR) {
                marker[i] = 2;
            }
        }
        
        if (marker[i] != a_bc_marker) {
            un1++;
            if (marker[i] != s2_bc_marker) un2++;
            if (marker[i] != s3_bc_marker) un3++;
            if (marker[i] != s2_bc_marker && marker[i] != s3_bc_marker) un23++;
        }

        bc_marker[i] = marker[i];  // Boundary marker
    }

    c1 = un1, c12 = un2, c13 = un3, c123 = un23;

    // Polarizations with Symmetry Bcs
    for( int k=0;k < nodes; k++ ) {
        if ( (!PEC_PEC) && polar==X_POLAR ) {
            // For X-polarization (PMC-PEC)
            bs2[k]         = bc_marker[k];  // Hx =0 along X-axis 
            bs2[nodes+k]   = -5;            // => Hy != 0 along X-axis @ bcm = 2 => PEC
            bs2[2*nodes+k] = bc_marker[k];  // Hz =0 along X-axis  @ bcm = 2

            bs3[k]         = bc_marker[k];  // => Hx = 0 along Y-axis@ bcm = 3
            bs3[nodes+k]   = -5;            // => Hy != 0 along Y-axis 
            bs3[2*nodes+k] = -5;            // Hz !=0 along Y-axis @ bcm = 3

            size = 2*c12 + c13;
        } else if ( (!PEC_PEC) && polar==Y_POLAR ) {
            // For Y-polarization (PEC-PMC)
            bs2[k]         = -5;            // Hx !=0 along X-axis 
            bs2[nodes+k]   = bc_marker[k];  // => Hy = 0 along X-axis @ bcm = 2 => PEC
            bs2[2*nodes+k] = -5;            // Hz !=0 along X-axis  @ bcm = 2
            
            bs3[k]         = -5;            // => Hx != 0 along Y-axis@ bcm = 3
            bs3[nodes+k]   = bc_marker[k];  // => Hy = 0 along Y-axis =>PMC  along Y-axis @ bcm = 3
            bs3[2*nodes+k] = bc_marker[k];  // Hz =0 along Y-axis @ bcm = 3

            size = c12 + 2*c13;
        } else if ( PEC_PEC && polar==X_POLAR) {
            // For X-polarization (PMC-PMC)
            bs2[k]         = bc_marker[k];  // Hx =0 along X-axis 
            bs2[nodes+k]   = -5;            // => Hy != 0 along X-axis @ bcm = 2 => PEC
            bs2[2*nodes+k] = bc_marker[k];  // Hz =0 along X-axis  @ bcm = 2

            bs3[k]         = -5;            // => Hx != 0 along Y-axis@ bcm = 3
            bs3[nodes+k]   = bc_marker[k];  // => Hy = 0 along Y-axis =>PMC  along Y-axis @ bcm = 3
            bs3[2*nodes+k] = bc_marker[k];  // Hz = 0 along Y-axis @ bcm = 3

            size= c12 + c13 + c123;
        } else if ( PEC_PEC && polar==Y_POLAR ) {
            // For Y-polarization (PEC-PEC)
            bs2[k]         = -5;            // Hx !=0 along X-axis 
            bs2[nodes+k]   = bc_marker[k];  // => Hy = 0 along X-axis @ bcm = 2 => PEC
            bs2[2*nodes+k] = -5;            // Hz !=0 along X-axis  @ bcm = 2

            bs3[k]         = bc_marker[k];  // => Hx = 0 along Y-axis@ bcm = 3 => PEC
            bs3[nodes+k]   = -5;            // => Hy != 0 along Y-axis 
            bs3[2*nodes+k] = -5;            // Hz !=0 along Y-axis @ bcm = 3

            size = c12 + c13 + c1 + ii;
        }

        // Artificial Bcs
        bs1[k]         = bc_marker[k];  // Hx =0
        bs1[nodes+k]   = bc_marker[k];  // Hy =0
        bs1[2*nodes+k] = bc_marker[k];  // Hz =0

        // Hz !=0 at (0, 0)
        if ( PEC_PEC && polar==Y_POLAR && ii!=0 && k==point0 )
            bs1[2*nodes+k]=0;
    }

    return 0;
}

int set_matrix(int je, int *p, double *x, double *y) {
    // distance of the elem from origin
    double centX=(x[1]+x[2]+x[0])/3.0;
    double centY=(y[1]+y[2]+y[0])/3.0;
    double d1=sqrt(centX*centX + centY*centY);

    //
    // RI distribution
    //
    double kp; // k_permitivity

    if  (elem_marker[je]==0)        kp=1.00;      // 0. Outside of PCF holy region
    else if (elem_marker[je]== 1)   kp=1/(n1*n1); // 1. Silica region of PCF
    else if (elem_marker[je]== 2)   kp=1/(n2*n2); // 2. Cladding Air-holes of PCF holey region
    else if (elem_marker[je]== 3)   kp=1/(n3*n3); // 3. core of PCF holey region

    double  a[3], b[3],c[3];

    for(int i=0;i<3;i++) {
        int j=(i+4)%3;
        int k=(i+5)%3;
        
        a[i]=x[j]*y[k]-x[k]*y[j];
        b[i]=y[j]-y[k];            // b = derivative  w.r.t  x   and
        c[i]=x[k]-x[j];            // c = derivative  w.r.t  y.
    }

    // AREA OF je-ELEMENT
    double area = fabs(0.5*(c[2]*b[1]-c[1]*b[2]));
    elem_area[je] = area;

    // i for row
    for(int i=0;i<3;i++) {
        for(int k=0;k<3;k++) {
            //AREA OF je-ELEMENT is included.
            NxNx[i][k]=b[i]*b[k]/(4*area);
            NyNy[i][k]=c[i]*c[k]/(4*area);
            NxNy[i][k]=b[i]*c[k]/(4*area);

            NNy[i][k] =c[k]/6.0;
            NNx[i][k] =b[k]/6.0;

            if (i == k) NN[i][k]=(area/6.0);
            else        NN[i][k]=(area/12.0);
        }
    }

    // Ref 1: Performance comparsion of FE approaches for EM waveguides
    // by  S.Selleri  and M. Zaboli , 1997
    double R[4][4];              // R[1]->R1, R[2]->R2, R[3]->R3
    double ps = 1/(0.5*(n1+n2)); // penalty term ; n2  is not a core index

    // i for row
    for(int i=0;i<3;i++) {
        for(int k=0;k<3;k++) {
            // Eq. 10 is calculated as under from REF (1) including penalty term ps.
            // -[S] - beta2* [T] =0 ....................(10) in Ref 1.
            //  Matrix S 
            R[1][1] =  -k02*NN[i][k]    + (kp*NyNy[i][k] + ps*NxNx[i][k]);  // kxx
            R[1][2] =  -kp *NxNy[i][k]  +  ps*NxNy[k][i];                   // kxy
            R[1][3] =   0;                                                  // kxz
            
            R[2][1] =  -kp *NxNy[k][i]  +  ps*NxNy[i][k];                   // kyx
            R[2][2] =  -k02*NN[i][k]    + (kp*NxNx[i][k] + ps*NyNy[i][k]);  // kyy
            R[2][3] =   0;                                                  // kyz
            
            R[3][1]=    kp *NNx[k][i]   +  ps*NNx[i][k];                    // kzx
            R[3][2] =   kp *NNy[k][i]   +  ps*NNy[i][k];                    // kzy
            R[3][3] =  -k02*NN[i][k]    + (kp*NyNy[i][k] + kp*NxNx[i][k]);  // kzz

            // Assembly of  [S]
            for ( int u=0; u < 3; u++ )
                for ( int l=0; l < 3; l++ )
                    S.set_element(p[i]+l*nodes, p[k]+u*nodes,
                                  S.find(p[i]+l*nodes, p[k]+u*nodes) + R[l+1][u+1]);
#ifdef DEBUG
            if (p[i]+2*nodes==22074) {
                for ( int u=0; u < 3; u++ )
                    printf("%d,%d,%f\t", p[i]+2*nodes, p[k]+u*nodes, S.find(p[i]+2*nodes, p[k]+u*nodes));
                printf("\n");
            }
#endif
            //  Matrix T
            R[1][1] =   kp*NN[i][k];
            R[1][2] =   0;
            R[1][3] =   kp*NNx[i][k] + ps* NNx[k][i];  // Row 1
            R[2][1] =   0;
            R[2][2] =   kp*NN[i][k];
            R[2][3] =   kp*NNy[i][k] + ps* NNy[k][i];  // Row 2
            R[3][1] =   0;
            R[3][2] =   0;
            R[3][3] =   ps*NN[i][k];                   // Row 3 

            // Assembly of  [T]
            for ( int u=0; u < 3; u++ )
                for ( int l=0; l < 3; l++ )
                    T.set_element(p[i]+l*nodes, p[k]+u*nodes,
                                  T.find(p[i]+l*nodes, p[k]+u*nodes) + R[l+1][u+1]);
        }
    }

    return 0;
}

int eigensystem() {
    n = size;    // Matrix Size of Eigen value problem
    int nnzA, nnzB;
    int *irowA = new int[NNZ_PER_COL*n];
    int *irowB = new int[NNZ_PER_COL*n];
    int *pcolA = new int[n+1];
    int *pcolB = new int[n+1];
    double *A = new double[NNZ_PER_COL*n];
    double *B = new double[NNZ_PER_COL*n];

    get_csc_matrix(S, nnzA, A, irowA, pcolA);
    get_csc_matrix(T, nnzB, B, irowB, pcolB);

    for ( int i=0; i < nnzA; i++ )
        A[i] = -1. * A[i];

    for ( int i=0; i < nnzB; i++ )
        B[i] = k02 * B[i];

    nev = nmodes;
    double sigma = n1*n1;
    EigValR = new double[nev*2];
    EigValI = new double[nev*2];
    EigVec = new double[nev*3*n];
#ifdef DEBUG
    FILE *fmat;
    fmat = fopen("matrixA.out", "w");
    for ( int k=0,i=0; i < n; i++ ) { 
        for ( int j=pcolA[i]; j < pcolA[i+1]; j++,k++)
            fprintf(fmat, "%d,%d,%f\t", i, irowA[k], A[k]);
        fprintf(fmat, "\n");
    }
    fclose(fmat);
    fmat = fopen("matrixB.out", "w");
    for ( int k=0,i=0; i < n; i++ ) { 
        for ( int j=pcolB[i]; j < pcolB[i+1]; j++,k++)
            fprintf(fmat, "%d,%d,%f\t", i, irowB[k], B[k]);
        fprintf(fmat, "\n");
    }
    fclose(fmat);
    exit(1);
#endif
    nconv = AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA, pcolA, 
                  nnzB, B, irowB, pcolB, sigma, nev);
#ifdef DEBUG
    printf("n %d\n", n);    
    printf("sigma %f\n", sigma);
    printf("NEV %d\n", nev);
    printf("nnzA %d\n", nnzA);
    printf("nnzB %d\n", nnzB);
    for ( int i=0; i < nev; i++)
        printf("%f + i %f\n", EigValR[i], EigValI[i]);
    for ( int i=0; i < 3*n*nev; i+=3)
        printf("%f\t%f\t%f\n", EigVec[i], EigVec[i+1], EigVec[i+2]);
    exit(1);
#endif
    delete [] irowA;
    delete [] irowB;
    delete [] pcolA;
    delete [] pcolB;
    delete [] A;
    delete [] B;

    return 0;
}

int get_csc_matrix(Matrix &mat, int &nnz, double *val, int *irow, int *pcol) {
    int index = 0;
    nnz = 0;

    pcol[0] = 0;
    for ( int col=0, colindex=0; col < 3*nodes; col++) {
        if (bs1[col]==a_bc_marker||bs2[col]==s2_bc_marker||bs3[col]==s3_bc_marker)
            continue;
        ++colindex;

        const int rowsize = mat.get_rowsize(col);
        int matindex = 0;
        int matrows[rowsize];
        double matvals[rowsize];

        mat.get_row_val(col, matrows, matvals);

        for ( int row=0, rowindex=-1; row < 3*nodes; row++ ) {
            if (bs1[row]==a_bc_marker||bs2[row]==s2_bc_marker||bs3[row]==s3_bc_marker)
                continue;
            ++rowindex;

            while ( matrows[matindex] < row )
                ++matindex;

            // end of row list
            if ( matindex >= rowsize )
                continue;

            if ( matrows[matindex] == row && matvals[matindex] != 0) {
                irow[index] = rowindex;
                val[index] = matvals[matindex];
                ++nnz;
                ++index;
            }
        }
        pcol[colindex] = index;
    }
    return 0;
}
