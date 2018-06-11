#ifndef __ENERGY_H
#define __ENERGY_H

#include "Graph.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

double Ecells(GLattice *Lattice);
void dEcells(GLattice *Lattice, double *dEdx, double *dEdy);

double Ebonds(GLattice *Lattice);
void dEbonds(GLattice *Lattice, double *dEdx, double *dEdy);

double E(const gsl_vector *x, void *params);
void  dE(const gsl_vector *x, void *params, gsl_vector *df);
void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

#endif