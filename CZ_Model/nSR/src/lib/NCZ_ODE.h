#ifndef NCZ_H
#define NCZ_H
#include <cstdlib>
#include <iostream>
// #include <cmath>
#include <math.h>

// #include "atrial_parameter.h"
#include "SingleCellParameter.hpp"
#include "Updated_CNZ.h"
void NCZ_ODE_Initialise(double *state, int celltype);
double NCZ_ODE(double dt, double Istim, SingleCellPara & het, double *state) ;


void NCZ_ODE_Update_Initialise(double *state, int celltype);
double NCZ_ODE_Update(double dt, double Istim, SingleCellPara & het, double *state) ;


void NCZ_ODE_New_Ito_Initialise(double *state, int celltype);
double NCZ_ODE_New_Ito(double dt, double Istim, SingleCellPara & het, double *state) ;


void NCZ_ODE_New_ItoItos_Initialise(double *state, int celltype);
double NCZ_ODE_New_ItoItos(double dt, double Istim, SingleCellPara & het, double *state);



// void NCZ_ODE_New_Itokv14_Initialise(double *state, int celltype, double ratio);

void NCZ_ODE_New_Itokv14_Initialise(double *state, int celltype, double ratio, int AF=0);
double NCZ_ODE_New_Itokv14(double dt, double Istim, SingleCellPara & het, double *state);
#endif