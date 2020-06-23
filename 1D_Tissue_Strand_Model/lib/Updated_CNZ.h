#ifndef CNZ_H
#define CNZ_H
#include <cstdlib>
#include <iostream>
// #include <cmath>
#include <math.h>

// #include "atrial_parameter.h"
#include "SingleCellParameter.hpp"

// #include "CNZ.tpp"
// template <class T>
double Updated_CNZ_ODE(double dt, double Istim, SingleCellPara & het, double *state) ;
double Updated_CNZ_ODE_Temp(double dt, double Istim, SingleCellPara & het, double *state) ;
// This is the actual model, that is probably called once per timestep.
// cell_params is a struct, defined in cell_params.h.  It contains (in v1) the
// current timestep, and scaling factors for the 13 main currents. While a cell
// doesn't have to use them, it does allow tweaking easily.
// state is an array of doubles, and contains the current state of the model.
// results is an array of doubles that can be used to carry transitory data into
// and out of the cell.
// double CNZ_ODE(double dt, double Istim, AtriaParameter *het, double *state) ;

void Updated_CNZ_ODE_Initialise(double *state, int celltype);


double INa_IKur_Drug_Updated_CNZ_ODE(double dt, double Istim, SingleCellPara & het, double *state);
void INa_IKur_Drug_Updated_CNZ_ODE_Initialise(double *state, int celltype);

int cnz_con_length(void);
// Returns the length of the state array


// double crn_con(cell_params *params, double *state, double *results)
// double
// The current state of the model is stored in state, and all but the membrane
// potential is updated during the timestep, using the euler method for
// concentrations and the rush-larsen method for the gating parameters.
// constants

double static const CRN_vcell = 20100.0; /* um3 */
double static const CRN_vi = 13668;
double static const CRN_vup = 1109.52/*0.0552*vcell*/;
double static const CRN_vrel = 96.48/*0.0048*vcell*/;
double static const T = 310; /* 37 Celcius */
double static const CRN_Tfac = 3;
double static const CRN_Csp = 1e+6; /* pF/cm2 */
double static const F = 96.4867; /* coul/mmol */
double static const R = 8.3143; /* J K-1 mol-1 */
double static const CRN_kb = 5.4; /* mM */
double static const CRN_nab = 140; /* mM */
double static const CRN_cab = 1.8; /* mM */
double static const CRN_nac = 140;
double static const CRN_cac = 1.8;
double static const CRN_kc = 5.4;
double static const CRN_gna = 7.8; /* nS/pF */
double static const CRN_gto = 0.1652; /* nS/pF */
double static const CRN_gkr = 0.8 * 0.029411765; /* nS/pF */
double static const CRN_gks = 0.8 * 0.12941176; /* nS/pF */
double static const CRN_gcaL = 0.1294 * 0.75; /* nS/pF */
double static const CRN_ErL = 65.0; /* mV */
double static const CRN_gk1 = 0.09; /* nS/pF */
double static const CRN_gbna = 0.8 * 0.0006744375; /* nS/pF */
double static const CRN_gbk = 0.00;
double static const CRN_gbca = 1.4 * 0.001131; /* nS/pF */
double static const CRN_inakbar = 1.4 * 0.59933874; /* pA/pF */
double static const CRN_kmnai = 10.0 * 0.95; /* mM */
double static const CRN_kmko = 1.5; /* mM */
double static const CRN_icapbar = 0.275; /* pA/pF */
double static const CRN_kmcap = 0.0005; /* mM */
double static const CRN_knacalr = 1600.0; /* pA/pF */
double static const CRN_kmnalr = 87.5; /* mM */
double static const CRN_kmcalr = 1.38; /* mM */
double static const CRN_ksatlr = 0.1;
double static const CRN_gammalr = 0.35;
double static const CRN_trpnbar = 0.070; /* mM */
double static const CRN_cmdnbar = 0.050; /* mM */
double static const CRN_csqnbar = 10; /* mM */
double static const CRN_kmcsqn = 0.8; /* mM */
double static const CRN_kmcmdn = 0.00238; /* mM */
double static const CRN_kmtrpn = 0.0005; /* mM */
double static const CRN_grelbar = 30.0; /* ms-1 */
double static const CRN_kmup = 0.00092; /* mM */
double static const CRN_iupbar = 0.005; /* mM/ms */
double static const CRN_caupmax = 15.0; /* mM */
double static const CRN_kupleak = 0.00033333336/*iupbar/caupmax*/; /* ms-1 */
double static const CRN_tautr = 180.0; /* ms */
double static const CRN_gf_na = 0.0944;////0.10/* 0.09*/; /* nS/pF */
double static const CRN_gf_k = 0.0752;
double static const CRN_gf = 0.025;//0.0752;
double static const CRN_Ef = -22.0; /* mV */

double static const gcaT = (0.15 * 0.22) / 17.62;
double static const EcaT = 45.0;

double static const MAL_gto = 0.75471 * 0.1962 * 1.0; // * 0.85; //nS/pF **** in manuscript it is given in nS, so has been /100 here
double static const MAL_gkur = 0.05874;
double static const CNZ_gkur = 0.006398 * 0.9 /** 1.1*/; //* 2.0*0.6;

double static const  GONG_gto = 0.103;
double static const  gf = 0.0385;


double static const Vss = 2 * 4.99232e-5;
double static const rjunct = 6.5;
double static const lcell = 122.051;
double static const dx = 1.625;


double static const Vnonjunct3 = 6 * 0.002531;
double static const VSR3 = 2 * 0.000057;
double static const VSR4 = 2 * 0.000080;
double static const Vcytosol = 0.008150;
// double static const Vnonjunct_Nai = 0.008100;

// Fibroblast stuff
double static const Cmf = 6.3; // pF
double static const vif = 0.00000137;//0.00000137;//0.00000000137; // um3  (from 0.00137 nL -> 0.00000137 nm3 -> 0.00000000137 um3
double static const naof = 130.0110; //mM
double static const kof = 5.3581; //mM
double static const Cm = 100.0;

double static const Aj_nj = 2492.324412;// = M_PI*rjunct*2*lcell*0.5; // atea between junct and non junct
double static const xj_nj = 0.822500;// = 0.02/2 + dx/2; // diffusion distance from center to center of junct to first njunct
// double static const xj_nj_Nai = 3.260000;// = 0.02/2 + 2*dx; // diffusion distance from center of junct to center of njunct (between 2nd and 3rd njunct)

double static const BCa = 24e-3;
double static const SLlow = 165;
double static const SLhigh = 13;

double static const KdBCa = 2.38e-3;
double static const KdSLlow = 1.1;
double static const KdSLhigh = 13e-3;
double static const CSQN =  6.7;
double static const KdCSQN = 0.8;
double static const BNa = 1.1319;
double static const KdBNa = 10;
double static const DCa = 780;// % µm2/s
double static const DCaSR = 44;//
double static const DCaBm = 25; //% µm2/s
double static const DNa = 0.12;
double static const SERCAKmf = 0.25e-3;
double static const SERCAKmr = 1.8;
double static const k4 = 7.5 ;// % pump rate
double static const k3 = 2.314815; //= k4 / SERCAKmr^2;
double static const k1 = 7500000.;// = 1000^2 * k4;
double static const k2 = 0.468750;// = k1 * SERCAKmf^2;
double static const cpumps = 40e-3;
double static const kSRleak = 6e-3;



#endif
