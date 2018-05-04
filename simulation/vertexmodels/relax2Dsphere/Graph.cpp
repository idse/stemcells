/*=================================================================
 *
 * GLattice data structure
 *
 *=================================================================*/

#include <math.h>
#include <stdlib.h>
#include "Graph.h"

// ------------------- Lattice constructor ----------------//

GLattice* Graph_Init(const mxArray *verts, const mxArray *bonds, const mxArray *cells,
                     const mxArray *scaleParam, const mxArray *T, const mxArray *A0,
                     const mxArray *cellPerim, const mxArray *l0, const mxArray *pT0,
                     const mxArray *kA0)
{
    double *vertsPtr;
    double *bondsPtr;
    double *CellBondsPtr;
    
    vertsPtr = mxGetPr(verts);
    bondsPtr = mxGetPr(bonds);
    
    GLattice *L = (struct GLattice *) malloc(sizeof(GLattice));
    
    /* initialize vertices
     * coordinates are stored in a vector
     * (x1, y1, z1, x2, y2, z2 ,..., x_n, y_n, z_n) */
    
    L->nverts = (int) mxGetM(verts);
    L->Verts = (GVert *) malloc(sizeof(GVert)*L->nverts);
    if (!L->Verts) mexErrMsgTxt("Can't allocate L->Verts!");

    for (int i=0; i < L->nverts; i++) {
        L->Verts[i].x = vertsPtr[i];
        L->Verts[i].y = vertsPtr[i+L->nverts];
        L->Verts[i].z = vertsPtr[i+L->nverts*2];
    }
    
    /* initialize bonds
     * convert from double to long
     * bonds is stored as
     * [v1 v2; v3 v4; ...] */
    
    L->nbonds = (int) mxGetM(bonds);
    L->Bonds = (GBond *) malloc(sizeof(GBond)*L->nbonds);
    if (!L->Bonds) mexErrMsgTxt("Can't allocate L->Bonds!");
    
    for (int i=0; i < L->nbonds; i++){
        L->Bonds[i].v1 = (int) bondsPtr[i]-1;
        L->Bonds[i].v2 = (int) bondsPtr[i+L->nbonds]-1;
        L->Bonds[i].c1 = (int) bondsPtr[i+L->nbonds*2]; /* no -1 here */
        L->Bonds[i].c2 = (int) bondsPtr[i+L->nbonds*3]; /* no -1 here */
    }
    
    /* initialize cells
     * cells is stored as {[b1 b2 b3 b4] [b2 b3 b4] ...}
     * bonds possess a direction and are ordered clockwise */
    
    L->ncells = (int) mxGetN(cells);
    L->Cells = (GCell *) malloc(sizeof(GCell)*L->ncells);
    if (!L->Cells) mexErrMsgTxt("Can't allocate L->Cells!");
    
    // initialize cell bonds and force
    for (int i=0; i < L->ncells; i++) {
        
        L->Cells[i].nbonds = (int) mxGetN(mxGetCell(cells, i));
        L->Cells[i].Bonds  = (int *)malloc(sizeof(int)*L->Cells[i].nbonds);
        if (!L->Cells[i].Bonds) mexErrMsgTxt("Can't allocate cell bonds!");
        
        // array for the force on each bond of cell
        L->Cells[i].force  = (double *)malloc(2*sizeof(double)*L->Cells[i].nbonds);

        CellBondsPtr = mxGetPr(mxGetCell(cells, i));
        for (int j=0; j < L->Cells[i].nbonds; j++){
            
            L->Cells[i].Bonds[j] = (int) CellBondsPtr[j]-1;

            // initialize force to zero
            L->Cells[i].force[j] = 0;  // Fx
            L->Cells[i].force[L->Cells[i].nbonds + j] = 0; // Fy
        }
    }

    /* Parameters */
    L->Paras = (double *) mxGetPr(scaleParam);
    L->A0 = (double *) mxGetPr(A0);
    L->kA0 = (double *) mxGetPr(kA0);
    L->cellPerim = (double *) mxGetPr(cellPerim);
    L->bondTension = (double *) mxGetPr(T);
    L->l0 = (double *) mxGetPr(l0);
    L->perimTension = (double *) mxGetPr(pT0);

    /* output */
    /* return the force as dEdx, dEdy, dEdz in one vector of length 3N  */
    L->tensionForce = (double *)malloc(3*sizeof(double)*L->nverts);
    L->pressureForce = (double *)malloc(3*sizeof(double)*L->nverts);
    L->perimeterForce = (double *)malloc(3*sizeof(double)*L->nverts);
    L->A = (double *)malloc(sizeof(double)*L->ncells);
    L->T = (double *)malloc(sizeof(double)*L->nbonds);

    // initialize force to zero
    for (int i=0; i < 3*L->nverts; i++) {
        L->tensionForce[i]      = 0;
        L->pressureForce[i]     = 0;
        L->perimeterForce[i]    = 0;
    }
    
    return L;
}

// ------------------- Vertex setter from GSL vector ----------------//

void GSLToVerts(const gsl_vector *x, GLattice *Lat)
{
    for (int i = 0; i<Lat->nverts; i++) {
        Lat->Verts[i].x = gsl_vector_get (x, 3*i);
        Lat->Verts[i].y = gsl_vector_get (x, 3*i+1);
        Lat->Verts[i].z = gsl_vector_get (x, 3*i+2);
    }
}

// ------------------- geometry of cells and bonds ----------------//

double CellArea(int icell, GLattice *L)
{
    double area, da;
    int ib;
    GVert v1, v2;
    
    GVert CM = CenterOfMass(L->Cells[icell], L);
    
    area = 0.;
    // mexPrintf("CM %f %f %f\n", CM.x, CM.y, CM.z);
    
    for (int i = 0; i <L->Cells[icell].nbonds; i++) {
        
        ib = L->Cells[icell].Bonds[i];
        v1 = L->Verts[ L->Bonds[ib].v1 ];
        v2 = L->Verts[ L->Bonds[ib].v2 ];
        
        // calculating areas relative to CM allows for roughly testing self intersection or extreme concavity
        da = (v1.x - CM.x)*(v2.y - CM.y) - (v2.x - CM.x)*(v1.y - CM.y);
        
        if (da < 0) {
            mexPrintf("-ve area contrib cell %d\n", icell);
        }
        area += da;
    }
    area *= 0.5;

    if (area < 0){
        mexPrintf("Negative area, bad vertex ordering?\n");//mexErrMsgTxt
        area = 0.0001;
    }

    return area;
}

double cellPerimeter(int icell, GLattice *L)
{
    double p = 0;
    int bi;
    
    for (int i = 0; i < L->Cells[icell].nbonds; i++){
        
        bi = L->Cells[icell].Bonds[i];
        p += BondLength(bi, L);
    }
    return p;
}

double BondLength(int ibond,struct GLattice *L)
{
    double length;
    int v1,v2;
    GVert *Verts = L->Verts;
    
    v1 = L->Bonds[ibond].v1;
    v2 = L->Bonds[ibond].v2;
    length = sqrt((Verts[v1].x-Verts[v2].x)*(Verts[v1].x-Verts[v2].x) +
                  (Verts[v1].y-Verts[v2].y)*(Verts[v1].y-Verts[v2].y));
    return length+0.0001;
}

GVert CenterOfMass(GCell C, struct GLattice *L){
    
    GVert CM, V;
    
    // initialize values to zero!
    CM.x = 0;   CM.y = 0;   CM.z = 0;
    
    for (int i = 0; i < C.nbonds; i++){
        
        V = L->Verts[L->Bonds[C.Bonds[i]].v1];
        CM.x += V.x;
        CM.y += V.y;
        CM.z += V.z;
    }
    CM.x = CM.x/C.nbonds;
    CM.y = CM.y/C.nbonds;
    CM.z = CM.z/C.nbonds;
    
    return CM;
}

// ------------------- Lattice deconstructor ----------------//

void Graph_Free(GLattice *L)
{
    for (int i=0; i < L->ncells; i++) {
        free(L->Cells[i].Bonds);
    }
    free(L->Cells);
    free(L->Bonds);
    free(L->Verts);
    free(L);
}
