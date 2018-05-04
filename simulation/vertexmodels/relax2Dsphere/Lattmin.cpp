/*=================================================================
 *
 * LattMin.C  .MEX file corresponding to LATTMIN.M
 *	          Minimizes energy of a lattice
 *
 * The calling syntax is:
 *
 *		[vertices] = latmin(g,par)
 *
 *=================================================================*/

#include <math.h>
#include "mex.h"
#include "energy.h"
#include "Graph.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#define DEBUG true

int MinimizeWithGSL(double *vertsout, GLattice *inLattice, int maxIter)
{
    size_t iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;    // gradient based (fdf) minimization algorithm
    gsl_multimin_fdfminimizer *s;               // minimizer state
    gsl_multimin_function_fdf M;                // function to minimize
    gsl_vector *xin, *xout;                     // initial and final vertex positions
    
    /* ------------------------------------------------------------------------------
     * initialize minimizer state, s, for algorithm T, model M, initial position xin
     * ------------------------------------------------------------------------------*/

    // copy vertex position from Lattice to xin
    xin  = gsl_vector_alloc(inLattice->nverts*3);
    for (int i = 0; i < inLattice->nverts; i++) {
        gsl_vector_set (xin, 3*i,   inLattice->Verts[i].x);
        gsl_vector_set (xin, 3*i+1, inLattice->Verts[i].y);
        gsl_vector_set (xin, 3*i+2, inLattice->Verts[i].z);
    }
    
    // the vertex model (dimension, initial condition, energy and energy gradient)
    M.n = xin->size;
    M.params = (void *) inLattice;
    M.f = &E;
    M.df = &dE;
    M.fdf = &my_fdf;
    
    // algorithm type; alternatives : _vector_bfgs, _conjugate_fr,_conjugate_pr
    // gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-with-Derivatives.html
    T = gsl_multimin_fdfminimizer_conjugate_pr;
    
    // minimizer of type T on space of dimension M.n
    s = gsl_multimin_fdfminimizer_alloc (T, M.n);
    
    // miminize M, starting at xin
    gsl_multimin_fdfminimizer_set(s, &M, xin, 2e-2, 2e-10); // 0.02, 2e-10); step size, tolerance
    
    /* -----------------------------------------------
     * update s using the iteration T
     * ----------------------------------------------- */
    
    status = GSL_CONTINUE;
    iter = 0;
    while (status == GSL_CONTINUE && iter < maxIter)
    {
        status = gsl_multimin_fdfminimizer_iterate (s);
        //mexPrintf("it %d\n", iter);
        
        // GSL_ENOPROG signifies that the minimizer is unable to improve
        if (status == GSL_ENOPROG){
            mexPrintf ("No progess: stop; iterations: %d \n", iter);
            break;
        }
        
        // tests the norm of the gradient g against the absolute tolerance epsab
        // return GSL_SUCCESS if |g| < epsabs
        status = gsl_multimin_test_gradient (s->gradient, 1e-5);
        if (status == GSL_SUCCESS) mexPrintf ("Minimum found after %d iterations\n", iter);
        
        iter++;
    };
    if (iter == maxIter) mexPrintf ("Reached max iterations: %d \n", iter);
    
    /* -----------------------------------------------
     * copy result to output variables and clean up
     * ----------------------------------------------- */
    
    xout  = gsl_vector_alloc(inLattice->nverts*3);
    gsl_vector_memcpy(xout,s->x);
    for (int i=0; i<inLattice->nverts; i++) {
        vertsout[i]                     = gsl_vector_get (xout, 3*i);
        vertsout[i+inLattice->nverts]   = gsl_vector_get (xout, 3*i+1);
        vertsout[i+inLattice->nverts*2] = gsl_vector_get (xout, 3*i+2);
    }
    
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free (xin);
    gsl_vector_free (xout);
    
    return  1;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{
    // define inputs
    # define VERTS      prhs[0]
    # define BONDS      prhs[1]
    # define CELLS      prhs[2]
    # define SCALEPARAM prhs[3]
    # define T0         prhs[4]
    # define A0         prhs[5]
    # define PERIM0     prhs[6]
    # define l0         prhs[7]
    # define PT0        prhs[8]
    # define kA0        prhs[9]

    int maxIter = (int)(*mxGetPr(prhs[10]));
    mexPrintf("max iterations %d\n", maxIter);
    
    // define outputs
    # define VERTSOUT       plhs[0]
    # define STATUS         plhs[1]
    # define ENERGY         plhs[2]
    # define TENSIONFORCE   plhs[3]
    # define PRESSUREFORCE  plhs[4]
    # define PERIMETERFORCE plhs[5]
    # define AREA           plhs[6]
    # define TENSION        plhs[7]
    # define CELLFORCE      plhs[8]
    
    // pointers to output data
    double *vertsout, *status, *energy, *tensionForce, *pressureForce, *perimeterForce, *area, *tension, *cf;
    int dims[2];
    
    // create main graph structure
    GLattice *Lattice = Graph_Init(VERTS, BONDS, CELLS, SCALEPARAM, T0, A0, PERIM0, l0, PT0, kA0);
    
    // initialize output variables
    STATUS = mxCreateDoubleMatrix(1, 1, mxREAL);
    status  = mxGetPr(STATUS);

    VERTSOUT = mxCreateDoubleMatrix((int) mxGetM(VERTS), 3, mxREAL);
    vertsout= mxGetPr(VERTSOUT);
    
    // Actual minimization using GSL
    status[0] = (double) MinimizeWithGSL(vertsout, Lattice, maxIter);
    
    // assign additional output variables from lattice
    ENERGY = mxCreateDoubleMatrix(1, 1, mxREAL);
    energy = mxGetPr(ENERGY);
    *energy = Lattice->energy;
    
    TENSIONFORCE = mxCreateDoubleMatrix((int) mxGetM(VERTS), 3, mxREAL);
    tensionForce = mxGetPr(TENSIONFORCE);
    
    PRESSUREFORCE = mxCreateDoubleMatrix((int) mxGetM(VERTS), 3, mxREAL);
    pressureForce = mxGetPr(PRESSUREFORCE);
    
    PERIMETERFORCE = mxCreateDoubleMatrix((int) mxGetM(VERTS), 3, mxREAL);
    perimeterForce = mxGetPr(PERIMETERFORCE);
    
    AREA = mxCreateDoubleMatrix(Lattice->ncells, 1, mxREAL);
    area = mxGetPr(AREA);
    
    TENSION = mxCreateDoubleMatrix(Lattice->nbonds, 1, mxREAL);
    tension = mxGetPr(TENSION);
    
    // return cell force
    dims[0] = Lattice->ncells;
    dims[1] = 1;
    CELLFORCE = mxCreateCellArray(1, dims);

    // there is probably a better way to do it than copying values
    for (int i = 0; i < Lattice->ncells; i++){
        mxArray *cellForceI = mxCreateDoubleMatrix(Lattice->Cells[i].nbonds, 2, mxREAL);
        cf = mxGetPr(cellForceI);
        mxSetCell(CELLFORCE, i, cellForceI);
        for (int j = 0; j < Lattice->Cells[i].nbonds; j++){
            cf[j] = Lattice->Cells[i].force[j];
            cf[Lattice->Cells[i].nbonds + j] = Lattice->Cells[i].force[Lattice->Cells[i].nbonds + j];
        }
    }

    for (int i = 0; i < 3*Lattice->nverts; i++){
        tensionForce[i]     = Lattice->tensionForce[i];
        pressureForce[i]    = Lattice->pressureForce[i];
        perimeterForce[i]   = Lattice->perimeterForce[i];
    }
    for (int i = 0; i < Lattice->ncells; i++){
        area[i]     = Lattice->A[i];
    }
    for (int i = 0; i < Lattice->nbonds; i++){
        tension[i]     = Lattice->T[i];
    }
    
    // Clean up
    Graph_Free(Lattice);
    
    return;  
}


