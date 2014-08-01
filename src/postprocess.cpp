#include "vfem.h"

extern int *bs1, *bs2, *bs3, *bc_marker, n;
extern double wavelength, n1, n2, n3, k0, k02;
extern POLAR polar;
extern int nmodes;
extern int elements, **id, *elem_marker, nodes;
extern double **coord, *elem_area;
extern int n, nconv, nev;
extern double *EigValR, *EigValI, *EigVec;

inline double sqr(double x) {
    return x*x;
}

// Marking the positions of unknown Ht-components for post processing
int post_process(std::string &output) {
    // noted: n = row2+1  also i.e the final matrix size
    int *nm = new int[3*nodes];
    for(int row=0, row2=-1; row < 3*nodes; row++) {
        if (bs1[row]==a_bc_marker||bs2[row]==s2_bc_marker||bs3[row]==s3_bc_marker)
            continue;
        row2++;
        nm[row2] = row;
    }

    FILE *f_out;
    
    f_out=fopen(output.c_str(), "w");

    double *Eff_index = new double[nconv];
    double *hx = new double[nodes];
    double *hy = new double[nodes];
    double *hz = new double[nodes];
    double **EV_Matrix = new double*[3*nconv];
    for (int i=0; i < 3*nconv; i++)
        EV_Matrix[i] = new double[nodes];

    // For eigen system serial
    fprintf(f_out,"Mode no: WaveLength (um), Eff-indices, Effective Mode area of quarter PCF (um2), ");
    fprintf(f_out,"Rel_sensitivity r %% for particular mode\n");

    for (int k=0; k < nconv ; k++) {
        Eff_index[k] = sqrt(EigValR[k]);
        printf("Eff_index[%2d] = %8.8f\n", k, Eff_index[k]);
        int beta1 = Eff_index[k]*k0;

        // n is the final matrix size ~ 3*nodes
        for (int i=0; i < n; i++) {
            if ( nm[i] < nodes ) {
                hx[nm[i]] = EigVec[k*n+i];
            } else if ( nm[i] < 2*nodes ) {
                hy[nm[i]-nodes] = EigVec[k*n+i];
            } else {
                hz[nm[i]-2*nodes] = EigVec[k*n+i];
            }
        }
        for (int i=0; i < nodes; i++) {
            EV_Matrix[k][i]         = hx[i]/beta1;
            EV_Matrix[k+nconv][i]   = hy[i]/beta1; // Plot Hy Field if polarization = 0 ELSE Plot Hx Field.
            EV_Matrix[k+2*nconv][i] = hz[i];
        }

        // IA=intensity in air, IS=Intensity in Silica, TP= Total insity in PCF holy region
        double IA, IS, TP, PF[nconv+1], Ia, Is;
        // Impedence of free space in   OHMS
        double Z0;
        // A_eff = Effective Mode area for PCF quarter structure  in  um^2
        double A_eff[nconv+1];
        double hx_avg, hy_avg, H, h2, H2, H4, h4;

        IA=0, IS=0, H2=0, H4=0, Z0=377.0;

        for (int je= 0; je < elements; je++) {
            hx_avg = (hx[id[je][0]]+hx[id[je][1]]+hx[id[je][2]])/3.;
            hy_avg = (hy[id[je][0]]+hy[id[je][1]]+hy[id[je][2]])/3.;
            H = sqrt(hx_avg*hx_avg + hy_avg*hy_avg); // Modulus of hx and hy.

            // 1. Effective Mode area
            // elem is inside holy PCF
            if (elem_marker[je] != 0) {
                h2=sqr(H)*elem_area[je];
                h4=sqr(H*H)*elem_area[je];
                H2+=h2;
                H4+=h4;
            }

            // 2. Relative Sensitivity r of particular mode
            // elem is in air hole = 2
            if ( elem_marker[je]  == 2 || (elem_marker[je] == 3 && n3 == n2) ) {
                Ia=elem_area[je]*sqr(H)*Z0/n2; // Intenisty in watt / um^2
                IA=IA+Ia;
            }
            // elem in silica region = 1
            if ( elem_marker[je]  == 1 || (elem_marker[je] == 3 && n3 == n1) )  {
                Is=elem_area[je]*sqr(H)*Z0/n1; // Intenisty in  watt / um^2
                IS=IS+Is;
            }
        }
        // Effective Mode area for PCF quarter structure in um^2
        A_eff[k]=sqr(H2)/H4;
        // Total power in Holy region of PCF in watt / um^2
        TP=IA+IS;
        // % Power fraction  FP => Rel_sensitivity r for particular mode in % =n2 * FP /(Eff_index[k])
        PF[k]=100.*IA/TP;
        fprintf(f_out,"%2d: %6.3f%15.9f%15.3f%12.2f (%)\n", 
                k, wavelength, Eff_index[k], A_eff[k], n2*PF[k]/Eff_index[k]);
    }

    
    // Printing Mode solution in array form (c3=Ist mode, c4= 2nd mode, ...etc) 
    fprintf(f_out,"Below columns (all).<< (left)Columns 2 & 3 are for x-coord and y-coord. ");
    fprintf(f_out,"Then Columns are in order of lowest to higher order modes.\n");
    for (int i=0; i < nodes; i++) {
        fprintf(f_out,"%5d[%d]:%20.15f%20.15f", i, bc_marker[i], coord[i][0], coord[i][1]);
        for (int k=0; k < 3*nconv; k++)
                fprintf(f_out, "%25.15f", EV_Matrix[k][i]); 
        fprintf(f_out,"\n");
    }
    fclose(f_out);

    delete [] hx;
    delete [] hy;
    delete [] hz;
    delete [] nm;
    delete [] Eff_index;
    for (int i=0; i < 3*nconv; i++)
        delete [] EV_Matrix[i];
    delete [] EV_Matrix;

    return 0;
}
