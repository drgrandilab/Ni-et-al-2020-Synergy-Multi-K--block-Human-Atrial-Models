#include "Lsoda_function_wrap.h"



void lsoda_generic_ODE(double t, double *Y, double * dY, void *usr_data) {


	SingleCellPara *Data = (SingleCellPara*) usr_data;


	// in the grandi model, Y[38] is the membrane potential
	// Y[38] = Data->V;
	Grandi_2011_ODE(t, Y, dY, Data );
	// Data->dV = dY[38];  // get the dV value save in the single cell data, and leave to solve using forward euler.
	dY[38] = 0.0;  // remember to set the derivative of vm in the single cell to be zero. 
}

