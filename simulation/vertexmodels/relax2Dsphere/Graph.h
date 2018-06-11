#ifndef __GRAPH_H
#define __GRAPH_H

#include <mex.h>
#include <gsl/gsl_multimin.h>

typedef struct GVert
{
    double x,y,z;       /* coordinates of vertex */
} GVert;

typedef struct GBond
{
    int v1,v2;          /* index of vertices */
    int c1,c2;          /* index of cells */
} GBond;

typedef struct GCell
{
    int nbonds;          /* number of cell */
    int *Bonds;          /* all bonds of a cell in clockwise order */
    double *force;
} GCell;

typedef struct GLattice
{
    // actual lattice
    
    int ncells;
    GCell *Cells;
    
    int nbonds;
    GBond *Bonds;
    
    int nverts;
    GVert *Verts;
    
    // parameters
    
    double *Paras;
    double *A0;
    double *kA0;
    double *cellPerim;
    double *bondTension;
    double *perimTension;
    double *l0;
    
    // quantities computed in lattmin and returned
    
    double energy;
    double *tensionForce;
    double *pressureForce;
    double *perimeterForce;
    double *A;
    double *T;
    
} GLattice;


GLattice* Graph_Init(const mxArray *verts, const mxArray *bonds, const mxArray *cells,
                     const mxArray *scaleParam, const mxArray *T, const mxArray *A0,
                     const mxArray *cellPerim, const mxArray *l0, const mxArray *pT0,
                     const mxArray *kA0);

void GSLToVerts(const gsl_vector *x, struct GLattice *Lat);
void Graph_Free(GLattice *Lattice);

double BondLength(int ibond, struct GLattice *Lattice);
double CellArea(int icell, struct GLattice *Lattice);
double cellPerimeter(int icell, GLattice *L);
GVert CenterOfMass(GCell C, struct GLattice *L);

#endif /* __GRAPH_H */