/*=================================================================
 *
 * energies for vertex model
 *
 *=================================================================*/

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include "Graph.h"

/* ----------------------- cell area energy and differential ----------------------------- */

// cell area energy
double EcellArea(GLattice *L){
    
    double *A0  = L->A0;
    double *kA0  = L->kA0;
    double c    = L->Paras[0];
    
    double energy = 0;

    for (int i=0; i < L->ncells; i++) {

        // area elasticity
        L->A[i] = CellArea(i,L);

        //energy += c*kA0[i]*(L->A[i] - A0[i])*(L->A[i]- A0[i])/2;
        energy += c*A0[i]/L->A[i];
    }
    
    //mexPrintf("cell energy %f\n", energy);
    
    return energy;
}

// cell area energy differential
void dEcellArea(GLattice *L, double *dEdx, double *dEdy){
    
    double *A0  = L->A0;
    double *kA0  = L->kA0;
    double c    = L->Paras[0];
    
    double A, prefactor;
    int vi, vpi, vmi, bimax;
    
    // gradient of area elasticity
    for (int i=0; i < L->ncells; i++) {
        
        A = CellArea(i, L);
        
        // loop over first element of each bond = loop over vertices around cell
        bimax = L->Cells[i].nbonds - 1;
        
        for(int j=0; j <= bimax; j++){

            vi = L->Bonds[ L->Cells[i].Bonds[j] ].v1;
            
            // cyclic previous vertex index
            if (j > 0)  { vmi = L->Bonds[ L->Cells[i].Bonds[j-1]    ].v1; }
            else        { vmi = L->Bonds[ L->Cells[i].Bonds[bimax]  ].v1; }
            
            vpi = L->Bonds[ L->Cells[i].Bonds[j] ].v2;
            
            // energy = (A-A0)^2/2
            //prefactor = c*kA0[i]*(A-A0[i])/2;
            
            // energy = 1/A
            prefactor = -c*A0[i]/(2*A*A);
            
            // update energy gradient
            dEdx[vi] += prefactor*(L->Verts[vpi].y - L->Verts[vmi].y);
            dEdy[vi] -= prefactor*(L->Verts[vpi].x - L->Verts[vmi].x);
            
            // update force on cell
            L->Cells[i].force[j] += -prefactor*(L->Verts[vpi].y - L->Verts[vmi].y);
            L->Cells[i].force[bimax + j + 1] += prefactor*(L->Verts[vpi].x - L->Verts[vmi].x);
        }
    }
            
    for (int i = 0; i < L->nverts; i++){
        // copy for output to matlab ; force = - gradient
        L->pressureForce[i]                 = -dEdx[i];
        L->pressureForce[i + L->nverts]     = -dEdy[i];
        L->pressureForce[i + 2*L->nverts]   = 0;
    }
}

/* ----------------------- cell perimeter energy and differential ----------------------------- */

// cell perimeter energy
double EcellPerim(GLattice *L){
    
    double *p0  = L->cellPerim;
    double *pT0 = L->perimTension;
    double c    = L->Paras[1];
    double k    = L->Paras[2];
    double alpha =  L->Paras[5];
    
    double energy = 0;
    double p;
    
    for (int ci = 0; ci < L->ncells; ci++) {
        
        p = cellPerimeter(ci, L);

        // perimeter elasticity
        energy += k*(p - p0[ci])*(p - p0[ci])/2;

        // perimeter string
        energy += c*pT0[ci]*pow(p,1+alpha);
    }

    return energy;
}

// cell perimeter energy differential
void dEcellPerim(GLattice *L, double *dEdx, double *dEdy){
    
    double *p0  = L->cellPerim;
    double *pT0 = L->perimTension;
    double c    = L->Paras[1];
    double k    = L->Paras[2];
    double alpha =  L->Paras[5];
    
    //mexPrintf("c %f k %f alpha %f", c,k,alpha);
    
    GVert *V    = L->Verts;
    
    double p, T, l;
    int cbi, bimax;
    int v1, v2;
    
    // reset force to zero
    for (int i=0; i < 3*L->nverts; i++) {
        L->perimeterForce[i] = 0;
    }
    
    // gradient of perimeter elasticity
    for (int ci = 0; ci < L->ncells; ci++) {
        
        p = cellPerimeter(ci, L);
        
        // loop over bonds
        bimax = L->Cells[ci].nbonds - 1;
        
        for(int bi=0; bi <= bimax; bi++){
            
            cbi = L->Cells[ci].Bonds[bi];

            // store the tension
            T = 0;
            T += k*(p - p0[ci]);
            T += c*pT0[ci]*(1+alpha)*pow(p,alpha);
            L->T[cbi] += T;
            
            // get the vertices involved in the bond
            v1 = L->Bonds[cbi].v1;
            v2 = L->Bonds[cbi].v2;
            
            // get the force
            l = BondLength(cbi, L); // the bond length
            L->perimeterForce[v1]               -= T*(V[v1].x-V[v2].x)/l;
            L->perimeterForce[v1 + L->nverts]   -= T*(V[v1].y-V[v2].y)/l;
            L->perimeterForce[v2]               += T*(V[v1].x-V[v2].x)/l;
            L->perimeterForce[v2 + L->nverts]   += T*(V[v1].y-V[v2].y)/l;
            
            // update force on cell
            L->Cells[ci].force[bi] -= T*(V[v1].x-V[v2].x)/l;
            if (bi < bimax) {
                L->Cells[ci].force[bi + 1] += T*(V[v1].x-V[v2].x)/l;
            }
            else{
                L->Cells[ci].force[0] += T*(V[v1].x-V[v2].x)/l;
            }
            L->Cells[ci].force[bi + bimax + 1] -= T*(V[v1].y-V[v2].y)/l;
            if (bi < bimax) {
                L->Cells[ci].force[bi + bimax + 2] += T*(V[v1].y-V[v2].y)/l;
            }
            else{
                L->Cells[ci].force[bimax + 1] += T*(V[v1].y-V[v2].y)/l;
            }
        }
    }
    
    // add this to the current value of the gradient
    for (int i = 0; i < L->nverts; i++){
        dEdx[i] += -L->perimeterForce[i];
        dEdy[i] += -L->perimeterForce[i + L->nverts];
    }
}

/* ----------------------- bond energy and differential ----------------------------- */

// bond energy
double Ebonds(GLattice *L){
    
//    double *restLengths = L->bondTension;
    double *T0   = L->bondTension;
    double *l0  = L->l0;
    double c    = L->Paras[3];
    double k    = L->Paras[4];
    double alpha =  L->Paras[5];
    
    double energy = 0;
    double bondenergy;
    double l;
    
    for (int i = 0; i < L->nbonds; i++){

        l = BondLength(i,L);
        bondenergy = 0;
        
        // bonds as Hookean springs
        bondenergy += k*(l-l0[i])*(l-l0[i])/2;
        
        // bonds as modified strings
        bondenergy += c*T0[i]*pow(l,1+alpha);
        
        // compensate for overcounting of internal bonds
        if      (L->Bonds[i].c2)  { bondenergy *= 0.5; }
        
        energy += bondenergy;
    }
    
    //mexPrintf("bond energy %f\n", energy);
    return energy;
}

// bond energy differential
void dEbonds(GLattice *L, double *dEdx, double *dEdy){
    
    double *T0   = L->bondTension;
    double *l0  = L->l0;
    GVert *V    = L->Verts;
    double c    = L->Paras[3];
    double k    = L->Paras[4];
    double alpha =  L->Paras[5];
    double l, T;
    int v1,v2;
    
    // reset force to zero
    for (int i=0; i < 3*L->nverts; i++) {
        L->tensionForce[i] = 0;
    }
    
    for (int i = 0; i < L->nbonds; i++){
        
        l = BondLength(i, L);
        
        T = 0;
        T += c*T0[i]*(1+alpha)*pow(l,alpha);    // tension of (generalized) strings
        T += k*(l - l0[i]);                     // tension of Hookean springs
        
        // compensate for overcounting of internal bonds
        if      (L->Bonds[i].c2)  { T *= 0.5; }
        
        // store for output
        L->T[i] += T;
        
        // get the vertices involved in the bond
        v1 = L->Bonds[i].v1;
        v2 = L->Bonds[i].v2;
        
        // add this to the current value of the gradient
        L->tensionForce[v1]             -= T*(V[v1].x-V[v2].x)/l;
        L->tensionForce[v1 + L->nverts] -= T*(V[v1].y-V[v2].y)/l;
        
        L->tensionForce[v2]             += T*(V[v1].x-V[v2].x)/l;
        L->tensionForce[v2 + L->nverts] += T*(V[v1].y-V[v2].y)/l;
    }
    
    // add -force to the current value of the gradient
    for (int i = 0; i < L->nverts; i++){
        
        dEdx[i] += -L->tensionForce[i];
        dEdy[i] += -L->tensionForce[i + L->nverts];
      }
}

/* ----------------------- total energy and differential ----------------------------- */

double E(const gsl_vector *x,void *params)
{
    // update vertex positions of Lattice
    GLattice *L = (GLattice *)params;
    GSLToVerts(x, L);

    // calculate energy
    double energy = 0;
    
    if (L->Paras[0] > 0)                    energy += EcellArea(L);
    if (L->Paras[1] > 0 || L->Paras[2] > 0) energy += EcellPerim(L);
    if (L->Paras[3] > 0 || L->Paras[4] > 0) energy += Ebonds(L);
    
    L->energy = energy;
    
    //mexPrintf("total energy %f\n", L->energy);
    
    return energy;
}

void dE(const gsl_vector *x, void *params, gsl_vector *df)
{
    GLattice *L = (GLattice *)params;
    double dEdx[L->nverts];
    double dEdy[L->nverts];
    double dEdz[L->nverts];
    
    // update vertex positions of Lattice
    GSLToVerts(x, L);
    
    // initialize gradient vector to zero
    for (int i=0; i < L->nverts; i++) {
        dEdx[i]=0;    dEdy[i]=0;    dEdz[i]=0;
    }

    // reset cell force to zero
    for (int i=0; i < L->ncells; i++) {
        for (int j=0; j< L->Cells[i].nbonds; j++){
            L->Cells[i].force[j] = 0;
            L->Cells[i].force[L->Cells[i].nbonds + j] = 0;
        }
    }
    
    // reset tension to zero
    for (int i=0; i < L->nbonds; i++) {
        L->T[i] = 0;
    }
    
    // update the gradient vector
    if (L->Paras[0] > 0)                    dEcellArea(L, dEdx, dEdy);
    if (L->Paras[1] > 0 || L->Paras[2] > 0) dEcellPerim(L, dEdx, dEdy);
    if (L->Paras[3] > 0 || L->Paras[4] > 0) dEbonds(L, dEdx, dEdy);

    // mexPrintf("call dEcellArea\n");
    
//    for (int i = 0; i < L->nverts; i++){
//        mexPrintf("total \t v %d x %.3f y %.3f \t dEdx %.2e dEdy %.2e\n", i, L->Verts[i].x, L->Verts[i].y, dEdx[i], dEdy[i]);
//    }
    
    // overwrite GSL output vector with updated gradient value
    for (int i = 0; i < L->nverts; i++) {
        gsl_vector_set(df, 3*i,   dEdx[i]);
        gsl_vector_set(df, 3*i+1, dEdy[i]);
        gsl_vector_set(df, 3*i+2, dEdz[i]);
    }
}

/* Compute both f and df together. */
void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f = E(x,params);
    dE(x,params,df);
}
