#ifndef __LATTMIN_H
#define __LATTMIN_H

int MinimizeWithGSL(double *vertsout, GLattice *inLattice, int maxIter);
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] );

#endif