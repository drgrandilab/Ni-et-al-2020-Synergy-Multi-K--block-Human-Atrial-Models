#ifndef RUDY_ISO_h
#define RUDY_ISO_h

// #include "CellTypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include "SingleCellParameter.hpp"


enum TypeCell { LVEPI, LVMCELL, LVENDO, RVEPI, RVMCELL, RVENDO, PFI, PFMB, PF};
enum IKsMutation { ORd, IKS_WT, HET, G269S};

static double const nao = 140.0; //extracellular sodium in mM
static double const cao = 1.8; //extracellular calcium in mM
static double const ko = 5.4; //extracellular potassium in mM

//buffer paramaters
static double const BSRmax = 0.047;
static double const KmBSR = 0.00087;
static double const BSLmax = 1.124;
static double const KmBSL = 0.0087;
static double const cmdnmax = 0.05;
static double const kmcmdn = 0.00238;
static double const trpnmax = 0.07;
static double const kmtrpn = 0.0005;
static double const csqnmax = 10.0;
static double const kmcsqn = 0.8;

//CaMK paramaters
static double const aCaMK = 0.05;
static double const bCaMK = 0.00068;
static double const CaMKo = 0.05;
static double const KmCaM = 0.0015;
static double const KmCaMK = 0.15;

//physical constants
static double const R = 8314.0;
static double const T = 310.0;
static double const F = 96485.0;


static double const L = 0.01;
static double const rad = 0.0011;
static double const vcell = 1000 * 3.14 * rad * rad * L;
static double const Ageo = 2 * 3.14 * rad * rad + 2 * 3.14 * rad * L;
static double const Acap = 2 * Ageo;
static double const vmyo = 0.68 * vcell;
static double const vmito = 0.26 * vcell;
static double const vsr = 0.06 * vcell;
static double const vnsr = 0.0552 * vcell;
static double const vjsr = 0.0048 * vcell;
static double const vss = 0.02 * vcell;
void initial_rudy_state(double *state, TypeCell celtype);
void INa_Drug_initial_rudy_state(double *state, TypeCell celtype);
double ORd_Model(double dt, double stims, double *state, TypeCell celltype, int ISO, IKsMutation mutation, float AB);
double INa_Drug_ORd_Model(double dt, double stims, double *state, TypeCell celltype, int ISO, IKsMutation mutation, float AB, SingleCellPara& para);

int ORd_con_length();
#endif