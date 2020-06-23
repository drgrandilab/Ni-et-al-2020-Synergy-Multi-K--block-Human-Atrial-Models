// Header file for the ICs and RHS functions for the Colman 2013 model

#ifndef CRN_1998_ODE_H
#define CRN_1998_ODE_H

// #include "atrial_parameter.h"
#include "SingleCellParameter.hpp"
// #include <cmath>
// #include "CRN.h"

void CRN_1998_ODE_Initialise(double *state, int celltype);
double CRN_1998_ODE(const double dt, const double Istim,  SingleCellPara &data, double *state);
double CRN_1998_ODE_New_INa(const double dt, const double Istim,  SingleCellPara &data, double *state);
#endif
