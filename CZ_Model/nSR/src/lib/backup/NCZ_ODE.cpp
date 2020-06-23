/*

 Implementation of IKur and INa state-dependent drug blockade. (using guarded-receptor model)

 This is the updatd CNZ model, update from the Colman et al. model 2013 Journal of Physiology Paper,

 Updates:  New IKur Model,
 Modifications to the conductances of several key ionic currents.

 Modifications to the Calcium dynamics model. change in the Refractory Parameters to make it faster, and thus
 remove the alternans seen at 400 ms of BCL.

 RMP: around -76 mV.

*/




#include "NCZ_ODE.h"
#include "ExplicitSolver.hpp"
#define CNZ_CON_STATE_LENGTH (40)






double NCZ_ODE(double dt, double Istim, SingleCellPara & het, double *state) {

	double V = state[0];
	double m = state[1];
	double h = state[2];
	double j = state[3];
	double d = state[4];
	double f = state[5];
	double xr = state[6];
	double xs = state[7];
	double nai = state[8];
	double cai = state[9];
	double ki = state[10];
	double fca = state[11];
	double itr = state[12];
	double its = state[13];
	double isusr = state[14];
	double isuss = state[15];
	double Cass = state[16];
	double CaSR1 = state[17];
	double CaSR2 = state[18];
	double SERCACa = state[19];
	double SERCACass = state[20];
	double RyRoss = state[21];
	double RyRcss = state[22];
	double RyRass = state[23];
	double RyRo3 = state[24];
	double RyRc3 = state[25];
	double RyRa3 = state[26];
	double dd = state[27];
	double ff = state[28];
	double CNZ_a = state[29];
	double CNZ_i = state[30];
	double If_y = state[31];
	double INa_BA = state[32];
	double INa_BI = state[33];
	double IKur_BO = state[34];
	double IKur_BC = state[35];
	double CNZ_is = state[36];



	double Ena, Ek, Eca, Ekf, Enaf;
	double alpha, beta, tau, inf, a, b;
	double INa, IKr, IKs, ICaL, IK1, IbK, IbCa;
	double IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT;
	double fnak, sigma;
	double naidot, kidot, caidot;


	// Reversal Potentials
	Ena = 26.71 * log(CRN_nac / nai);
	Ek = 26.71 * log(CRN_kc / ki);
	Eca = 13.35 * log(CRN_cac / cai);



	//INa
	INa = het.GNa * Cm * CRN_gna * (1 - INa_BA - INa_BI) * m * m * m * h * j * (V - Ena); // drug blocked INa

	double dBAdt = het.drug_INa_Ka * het.drug_INa_concen * m * m * m * h * j * ( 1 - INa_BA - INa_BI) - het.drug_INa_La * INa_BA ; // * exp(-z*V*F/R/T)
	double dBIdt = het.drug_INa_Ki * het.drug_INa_concen * (1 - h) * (1 - INa_BA - INa_BI) - het.drug_INa_Li * INa_BI; // * exp(-z*V*F/R/T)


	INa_BA += dt * dBAdt;   // blockade to Open gate.
	INa_BI += dt * dBIdt;	// Blockade to closed states.


	// m gate
	alpha =  0.32 * ((V - het.INa_Act_Shift) + 47.13) / (1 - exp(-0.1 * ((V - het.INa_Act_Shift) + 47.13)));
	if (fabs(V + 47.13) < 1e-10) alpha = 3.2;
	beta = 0.08 * exp(-(V - het.INa_Act_Shift) / 11.0);

	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	m = inf + (m - inf) * exp(-dt / tau); // steady state approx due to fast tau
	// m = inf + (m - inf) * exp(-dt / tau); // steady state approx due to fast tau

	// j gate
	if (V >= -40.0)
	{
		alpha  = 0.0;
		beta = 0.3 * exp(-2.535e-7 * V) / (1 + exp(-0.1 * (V + 32)));
	}
	else
	{
		alpha = (-1.2714e5 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) * (V + 37.78) / (1 + exp(0.311 * (V + 79.23)));
		beta = 0.1212 * exp(-0.01052 * V) / (1 + exp(-0.1378 * (V + 40.14)));
	}
	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	j = inf + (j - inf) * exp(-dt / tau);
	//j = table_row[INA_J_INF] + (j - table_row[INA_J_INF])*table_row[INA_J_TAU];

	// h gate
	if (V >= -40.0)
	{
		alpha  = 0.0;
		beta = 1 / (0.13 * (1 + exp((V + 10.66) / -11.1)));
	}
	else
	{
		alpha = 0.135 * exp((V + 80) / -6.8);
		beta = 3.56 * exp(0.079 * V) + 3.1e5 * exp(0.35 * V);
	}
	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	h = inf + (h - inf) * exp(-dt / tau);
	//h = table_row[INA_H_INF] + (h - table_row[INA_H_INF])*table_row[INA_H_TAU];

	// IKr
	IKr = 1.25 * het.GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4)); // 90% haibo
	// IKr = 0.9 * het.GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4)); // 90% haibo

	// xr gate
	a = 0.0003 * (V + 14.1) / (1 - exp((V + 14.1) / -5));
	b = 0.000073898 * (V - 3.3328) / (exp((V - 3.3328) / 5.1237) - 1);
	if (fabs(V + 14.1) < 1e-10) a = 0.0015; /* denominator = 0 */
	if (fabs(V - 3.3328) < 1e-10) b = 3.7836118e-4; /* denominator = 0 */
	tau = 1 / (a + b);
	inf = 1 / (1 + exp((V + het.IKr_ac_shift + 14.1) / (-6.5 * het.IKr_grad)));
	xr = inf + (xr - inf) * exp(-dt / tau);

	// IKs
	IKs = 1.25 * het.GKs * Cm * CRN_gks * xs * xs * (V - Ek); // 1.1 times larger. haibo

	double pNaK = 0.01833;
	double eks = 26.71 * log((CRN_kc + pNaK * CRN_nac) / (ki + pNaK * nai));
	// IKs = het.GKs * Cm * CRN_gks * xs * xs * (V - eks);  // use grandi style of reversal potential
	// xs gate
	a = 0.00004 * (V - 19.9) / (1 - exp((V - 19.9) / -17));
	b = 0.000035 * (V - 19.9) / (exp((V - 19.9) / 9) - 1);
	if (fabs(V - 19.9) < 1e-10) // denominator = 0
	{
		a = 0.00068;
		b = 0.000315;
	}
	tau = 0.5 / (a + b); // note lagrer taus may be more accurate
	inf = sqrt(1 / (1 + exp((V - 19.9 - het.IKs_shift) / (-12.7 * het.IKs_grad))));
	xs = inf + (xs - inf) * exp(-dt / tau);
	//xs = table_row[IKS_INF] + (xs - table_row[IKS_INF])*table_row[IKS_TAU];
	//ICaL
	// ICaL = 2.125 * het.GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL);
	// ICaL = 1.75 * het.GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL);   // haibo 1.75
	ICaL = 1.0 * 1.333333 * het.GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL); // haibo 1.75

	// fca gate
	inf = 1 / (1 + (Cass / 0.00035));
	tau = 2.0;
	fca = inf + (fca - inf) * exp(-dt / tau);

	// d gate
	a = 1 / (1 + exp((V + 10) / -6.24));
	tau = a * (1 - exp((V + 10) / -6.24)) / (0.035 * (V + 10));
	if (fabs(V + 10) < 1e-10) tau = a * 4.579;
	inf = 1 / (1 + exp((V + 10) / -8));
	d = inf + (d - inf) * exp(-dt / tau);

	// f gate
	inf = exp(-(V + 28) / 6.9) / (1 + exp(-(V + 28) / 6.9));
	tau = 1.5 * 2 * 3 / (0.0197 * exp(-0.0337 * 0.0337 * (V + 10) * (V + 10)) + 0.02);
	f = inf + (f - inf) * exp(-dt / tau);

	//Ito
	Ito = 1.08 * 1.0 * het.Gto * Cm * MAL_gto * itr * its * (V - Ek);

	//r gate
	// inf = 1.0 / (1.0 + exp(((V - het.Ito_Act_Shift) - 1.0) / -11.0));
	inf = 1.0 / (1.0 + exp(((V - het.Ito_Act_Shift) - 1.0) / -11.0));
	// tau = (0.0035 * exp(-((V - het.Ito_Act_Shift) / 30.0) * 2) + 0.0015);
	tau = 3.5 * ( exp(- pow((V / 30.0), 2.0))) + 1.5;
	// tau = 10 * ( exp(- pow(((V+40) / 45.0), 2.0))) + 0.79;  // change to 10 works...
	itr = inf + (itr - inf) * exp(-(dt) / tau);
	//itr = table_row[ITO_ITR_INF] + (itr - table_row[ITO_ITR_INF])*table_row[ITO_ITR_TAU];

	// s gate
	inf = 1.0 / (1.0 + exp((V + 40.5) / 11.5));
	tau = 1000.0 * (0.025635 * exp (-((V + 52.45) / 15.8827) * ((V + 52.45) / 15.8827)) + 0.01414);
	// tau = 25.635 * exp(-(pow((V + 52.45) / 15.8827, 2.0))) + 14.14;
	// tauytof = 85.0 * exp(-pow(Vm + 52.45, 2.0) / 15.8) + 24.14; // check haibo wrong version
	// tauytof = 25.635 * exp(-(pow((Vm + 52.45) / 15.8827, 2.0))) + 24.14; //%14.14  from matlab version from Grandi model and Malerckar et al.
	its = inf + (its - inf) * exp(-(dt) / tau);
	//its = table_row[ITO_ITS_INF] + (its - table_row[ITO_ITS_INF])*table_row[ITO_ITS_TAU];

	// 	//IKur
	// 	double IKur_K0 = -8.26597;
	// double IKur_c = 4.5128;
	// double IKur_x0 = 1.899769;
	// double IKur_y0 = 20.5232;


	IKur = 1.05 * 1.0 / 0.9 * het.GKur * Cm * CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597)))) * (1 - IKur_BO - IKur_BC) * CNZ_a * (0.65 * CNZ_i + 0.35 * CNZ_is) * (V - Ek);
	// IKur = het.GKur * Cm * CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597)))) * 0.5 * CNZ_a * CNZ_i * (V - Ek);

	double d_IKur_BO = het.drug_IKur_KO * exp(-het.drug_IKur_ZKO * V * F / R / T) * het.drug_IKur_concen * CNZ_i * CNZ_a * (1 - IKur_BO - IKur_BC) - het.drug_IKur_LO * IKur_BO * exp(-het.drug_IKur_ZLO * V * F / R / T);
	double d_IKur_BC = het.drug_IKur_KC * exp(-het.drug_IKur_ZKC * V * F / R / T) * het.drug_IKur_concen * CNZ_i * (1 - CNZ_a) * (1 - IKur_BO - IKur_BC) - het.drug_IKur_LC * IKur_BC * exp(-het.drug_IKur_ZLC * V * F / R / T);
	IKur_BO += dt * d_IKur_BO;
	IKur_BC += dt * d_IKur_BC;
	// IKur = GKur * IKur;

	/*inf = 1 / (1 + exp(-(V + IKur_ac1_shift + 6) / (8.6 * IKur_ac1_grad)));
	tau = (0.009 / (1.0 + exp((V + 5) / 12)) + 0.0005);
	isusr = inf + (isusr - inf) * exp(-(dt / 1000) / tau);

	inf = 1 / (1 + exp((V + IKur_inac_shift + 7.5) / (10 * IKur_inac_grad)));
	tau = (0.59 / (1 + exp((V + 60.0) / 10)) + 3.05);
	isuss = inf + (isuss - inf) * exp(-(dt / 1000) / tau);*/

	double K_Q10 = 3.0;//3.5308257834747638;
	//CNZ_a
	// inf = ((het.IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + het.IKur_ac1_shift) / (-5.75418 * het.IKur_ac1_grad))) ) * ((het.IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + het.IKur_ac2_shift) / (-11.51037561 * het.IKur_ac2_grad)))) + het.IKur_ac_add;


	inf = 1.0 / (1 + exp(-(V  + 5.52) / (8.6))); // change to simple version, haibo.  removed IKur mutations here.
	tau = (45.67 / (1.0 + exp((V + 9.23) / 10.53)) + 4.27) * ( 0.55 / (1.0 + exp((V + 60.87) / -22.88)) + 0.03); // incorporate deactivation Thu 20 Jul 2017 23:24:01 BST Haibo
	tau = tau / K_Q10;


	CNZ_a = inf + (CNZ_a - inf) * exp(-(dt) / tau);

	// CNZ_i
	inf = 1.0 / (1.0 + exp( (V + 7.5) / 6.67)); // Thu 20 Jul 2017 23:24:21 BST
	tau = 1300.0 / (exp((V - 133.0) / 97.0) * 5.04 + 0.2 * exp((V + 10.9233071) / (-12.0))) + 370.0; // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998


	CNZ_i = inf + (CNZ_i - inf) * exp(-(dt) / tau);
	tau = 1400.0 / (exp((V + 0.1) / 8.86) + 1.0) + 478.0 / (exp((V - 41.0) / -3.0) + 1.0) + 6500.0;// // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998

	CNZ_is = inf + (CNZ_is - inf) * exp(-(dt) / tau);

	//If

	double If_Na, If_K;
	If_Na = 0.0;//1.5*0.3833 * het.Gf * Cm * CRN_gf * If_y * (V - Ena);
	If_K =  0.0;//1.5* 0.6167 * het.Gf * Cm * CRN_gf * If_y * (V - Ek);

	// If = 1.0  * het.Gf * Cm * CRN_gf * If_y * (V - (-22));
	If = 1.0  * het.Gf * Cm * CRN_gf * If_y * (V - (-22));

	// If y gate
	inf = 1 / (1 + exp((V + 90.95 + het.If_vshift) / (10.1 * het.If_grad))   );
	tau = 1 / (1.2783E-04 * exp(-V / 9.2424) + 121.6092 * exp(V / 9.2424)   );
	If_y = inf + (If_y - inf) * exp(-dt / tau);

	//ICaT
	ICaT = 0.0; //  * het.GCaT * Cm * gcaT * ff * dd * (V - EcaT);

	// dd gate
	/*a = 1.0680 * exp((V + 26.3) / 30.0);
	b  = 1.0680 * exp(-(V + 26.3) / 30.0);
	tau = 1.0 / (a + b);
	inf = 1.0 / (1.0 + exp(-(V + 37.0) / 6.8));
	dd = inf + (dd - inf) * exp(-dt / tau);

	// ff gate
	a = 0.0153 * exp(-(V + 71.0) / 83.3);
	b  = 0.0150 * exp((V + 71.0) / 15.38);
	tau = 1.0 / (a + b);
	inf = 1.0 / (1.0 + exp((V + 71.0) / 9.0));
	ff = inf + (ff - inf) * exp(-dt / tau);*/

	//IK1
	IK1 =  het.GK1 * Cm * CRN_gk1 * (V - Ek + het.IK1_v_shift) / (1 + exp(0.07 * (V - (Ek + 6.94) + het.IK1_v_shift)));
	// IK1 =  het.GK1 * Cm * CRN_gk1 * (V - Ek + het.IK1_v_shift) / (1 + exp(0.07 * (V +80 + het.IK1_v_shift)));
	//IK1 = GK1 * Cm * CRN_gk1 * (V-Ek+het.IK1_v_shift) * (het.IK1_type_default * (table_row[IK1_FAC]) + het.IK1_type_AS *(table_row[IK1_FAC_AS])  );
	// printf("%f\n",  Ek);
	//Iab
	// Iab = 0.0;//0.* Cm * 0.0003879 * (V + 69.6) / (1 - 0.8377 * exp((V + 49.06) / 1056));

	//Ibk
	IbK =  0.0;//Cm * CRN_gbk * (V - Ek);

	//IbCa
	IbCa = 0.9 * Cm * CRN_gbca * (V - Eca); // increased background IbCa, haibo
	// IbNa  1.7
	IbNa = 1.25 * 1.0 * Cm * CRN_gbna * (V - Ena); // increase. haibo

	//ICap  // 1.26 // 1.15
	ICap = 1.0 * het.GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);  // change, haibo
	// ICap = 1.1 * het.GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);  // change, haibo



	//INaCa, times 1.4, Haibo
	INaCa = 1.1 * het.GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
	        (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
	        (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
	         exp(V * (CRN_gammalr - 1) * F / (R * T)));


	// Version 3.0`
	/*INaCa = 1.6 * het.GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
	        (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
	        (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
	         exp(V * (CRN_gammalr - 1) * F / (R * T)));*/


	//INaK
	sigma = (exp(CRN_nac / 67.3) - 1) / 7.0;
	fnak = 1.0 / (1.0 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0365 * sigma * exp(-V * F / (R * T)));
	// INaK = 1.3 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(/*CRN_kmnai*/11/ nai, 4));  // increase, Haibo.

	// verison 3.0  // 1.28
	INaK = 1.28 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(/*CRN_kmnai*/10.1/ nai, 4)); // increase, Haibo.

	// INaK = 1.28 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(CRN_kmnai / nai, 1.5));

	// concentrations
	naidot = (-3 * INaK - 3 * INaCa - IbNa - INa -  If_Na) / (F * CRN_vi);
	kidot = (2 * INaK - IK1 - Ito - IKur - IKr - IKs - IbK - If_K) / (F * CRN_vi);
	nai = nai + dt * naidot;
	// ki = ki + dt * kidot; // K+ is drifting

	//Calcium Handling
	double betass;
	betass = pow(( 1 + SLlow * KdSLlow / pow(((Cass) + KdSLlow), 2) + SLhigh * KdSLhigh / pow(((Cass) + KdSLhigh), 2) + BCa * KdBCa / pow(((Cass) + KdBCa), 2)  ), (-1));

	double betai, gammai;
	betai = pow(( 1 + BCa * KdBCa / pow((cai + KdBCa), 2)  ), (-1));
	gammai = BCa * KdBCa / pow((cai + KdBCa), 2);

	double betaSR1, betaSR2;
	betaSR1 = pow( ( 1 + CSQN * KdCSQN / pow((CaSR1 + KdCSQN), 2) ), (-1));
	betaSR2 = pow( ( 1 + CSQN * KdCSQN / pow((CaSR2 + KdCSQN), 2) ), (-1));


	double Jj_nj;
	Jj_nj = 2.5 * DCa * Aj_nj / xj_nj * ((Cass) - cai) * 1e-6;  // increase, Haibo

	double J_SERCASR, J_bulkSERCA;
	double J_SERCASRss, J_bulkSERCAss;

	J_SERCASR =    het.SERCA_Scale * 0.75 * (-k3 * pow(CaSR1, 2) * (cpumps - SERCACa) + k4 * (SERCACa)) * Vnonjunct3 * 2;
	J_bulkSERCA =   het.SERCA_Scale * 0.75 * (k1 * pow(cai, 2) * (cpumps - (SERCACa)) - k2 * (SERCACa)) * Vnonjunct3 * 2; // increase, haibo

	J_SERCASRss =   het.SERCA_Scale * 0.75 * (-k3 * pow(CaSR2, 2) * (cpumps - (SERCACass)) + k4 * (SERCACass)) * Vss * 2;
	J_bulkSERCAss =    het.SERCA_Scale * 0.75 * (k1 * pow((Cass), 2) * (cpumps - (SERCACass)) - k2 * (SERCACass)) * Vss * 2; // increase, haibo



	// change parameters, to remove alternans.
	double RyRtauadapt = het.RyRTauScale * 0.25;
	double RyRtauactss = 5e-3;
	double RyRtauinactss =  het.RyRTauScale * 15e-3;
	double RyRtauact =  5e-3;//0.6 * 8.75e-3;
	double RyRtauinact = het.RyRTauScale * 2 * 15e-3; //0.6 * 87.5e-3;

	// original parameters.
	/*double RyRtauadapt = 1;

	double RyRtauactss = 5e-3;
	double RyRtauinactss = 15e-3;
	double  RyRtauact = 18.75e-3;
	double RyRtauinact = 87.5e-3;*/

	double nuss = 625 * Vss;
	double RyRSRCass = (1 - 1 / (1 +  exp((CaSR2 - 0.3) / 0.1)));
	double RyRainfss = 0.505 - 0.427 / (1 + exp((( Cass + (het.fRyR * Cass) ) * 1000 - 0.29) / 0.082));
	double RyRoinfss = (1 - 1 / (1 +  exp(((Cass + (het.fRyR * Cass)   ) * 1000 - ((RyRass) + 0.22)) / 0.03)));
	double RyRcinfss = (1 / (1 + exp(((Cass + (het.fRyR * Cass )) * 1000 - ((RyRass) + 0.02)) / 0.01)));
	double Jrelss = nuss * ( (RyRoss) ) * (RyRcss) * RyRSRCass * ( CaSR2 -  (Cass) );

	double nu3 = Vnonjunct3;
	double RyRSRCa3 = (1 - 1 / (1 +  exp((CaSR1 - 0.30) / 0.1)));
	double RyRainf3 =  0.505 - 0.427 / (1 + exp(( (cai * 2 + ( het.fRyR * cai * 2)  ) * 1000 - 0.29) / 0.082));     // increased sensitiviy haibo
	double RyRoinf3 = (1 - 1 / (1 +  exp(( (cai * 2 + ( het.fRyR * cai * 2) ) * 1000 - ((RyRa3) + 0.22)) / 0.03)));   // increased sensitiviy haibo
	double RyRcinf3 = (1 / (1 +  exp(( (cai * 2 + (het.fRyR * cai * 2 ) ) * 1000 - ((RyRa3) + 0.02)) / 0.01)));   // increased sensitiviy haibo
	double Jrel3 = nu3 * ( (RyRo3) ) * (RyRc3) * RyRSRCa3 * ( CaSR1 -  cai );   // increased sensitiviy haibo

	Jrelss = het.fIRel * Jrelss;
	Jrel3 = het.fIRel * Jrel3;

	double JSRCaleak3 = 0.5 * het.GSR_leak * kSRleak * ( CaSR1 - cai ) * Vnonjunct3; // decreased leakage.  haibo
	double JSRCaleakss = 0.5 * het.GSR_leak * kSRleak * ( CaSR2 - (Cass) ) * Vss; // decreased leakage.  haibo

	double JCa, JCass;
	JCa = -/*het.BULK_CONST * */J_bulkSERCA + JSRCaleak3 + Jrel3 + Jj_nj;   // no bulk const haibo
	JCass = -Jj_nj + JSRCaleakss - /*het.BULK_CONST * */J_bulkSERCAss + Jrelss;  // remove bulk const haibo

	double JSRCa1, JSRCa2;
	JSRCa1 = J_SERCASR - JSRCaleak3 - Jrel3;
	JSRCa2 = J_SERCASRss - JSRCaleakss - Jrelss;

	double dy;

	dy = /*0.5 * */(-J_SERCASR + J_bulkSERCA) / Vnonjunct3;   // 0.5* does not make sense to me, remove it. haibo
	SERCACa = Foward_Euler(dt / 1000, dy, SERCACa);
	dy = /*0.5 * */(-J_SERCASRss + J_bulkSERCAss) / Vss;   // 0.5* does not make sense to me, remove it. haibo
	SERCACass = Foward_Euler(dt / 1000, dy, SERCACass);

	RyRoss = Euler_inf(dt / 1000, RyRoss, RyRoinfss, RyRtauactss);
	RyRcss = Euler_inf(dt / 1000, RyRcss, RyRcinfss, RyRtauinactss);
	RyRass = Euler_inf(dt / 1000, RyRass, RyRainfss, RyRtauadapt);
	RyRo3 = Euler_inf(dt / 1000, RyRo3, RyRoinf3, RyRtauact);
	RyRc3 = Euler_inf(dt / 1000, RyRc3, RyRcinf3, RyRtauinact);
	RyRa3 = Euler_inf(dt / 1000, RyRa3, RyRainf3, RyRtauadapt);

	dy =  betass * ( JCass / Vss + ((-( /*het.RyR * */ICaL) - IbCa - ICap - ICaT + 2 * INaCa)) / (2 * Vss * 1000 * F) );  // this parameter does not make sense to me, het.RyR , haibo
	Cass = Foward_Euler(dt / 1000, dy, Cass);

	dy = JCa / Vnonjunct3 * betai;
	cai = Foward_Euler(dt / 1000, dy, cai);

	dy =  betaSR1 * (DCaSR) * ( (CaSR2 - 2 * CaSR1 + CaSR1) / (dx * dx) + (CaSR1 - CaSR2) / (2 * 3 * (dx * dx)) ) + JSRCa1 / VSR3 * betaSR1;

	CaSR1 = Foward_Euler(dt, dy, CaSR1);

	dy = betaSR2 * (DCaSR) * ( (CaSR2 - 2 * CaSR2 + CaSR1) / (dx * dx) + (CaSR2 - CaSR1) / (2 * 4 * (dx * dx)) ) + JSRCa2 / VSR4 * betaSR2;

	CaSR2 = Foward_Euler(dt, dy, CaSR2);

	// return state variables
	state[1] = m;
	state[2] = h;
	state[3] = j;
	state[4] = d;
	state[5] = f;
	state[6] = xr;
	state[7] = xs;
	state[8] = nai;
	state[9] = cai;
	state[10] = ki;
	state[11] = fca;
	state[12] = itr;
	state[13] = its;
	state[14] = isusr;
	state[15] = isuss;
	state[16] = Cass;
	state[17] = CaSR1;
	state[18] = CaSR2;
	state[19] = SERCACa;
	state[20] = SERCACass;
	state[21] = RyRoss;
	state[22] = RyRcss;
	state[23] = RyRass;
	state[24] = RyRo3;
	state[25] = RyRc3;
	state[26] = RyRa3;
	state[27] = dd;
	state[28] = ff;
	state[29] = CNZ_a;
	state[30] = CNZ_i;
	state[31] = If_y;
	state[32] = INa_BA;
	state[33] = INa_BI;
	state[34] = IKur_BO;
	state[35] = IKur_BC;
	state[36] =  CNZ_is;


	het.INa = INa / Cm;
	het.IKur = IKur / Cm;
	het.IKr = IKr / Cm;
	het.ICaL = ICaL / Cm;
	het.Ito = Ito / Cm;
	het.IKs = IKs / Cm;
	het.INaK = INaK / Cm;
	het.IK1 = IK1 / Cm;
	het.INaCa = INaCa / Cm;

	return (Istim + (INa + Ito + IKur + IKr + IKs + ICaL + IK1 + IbK + IbNa + IbCa + INaCa + INaK + ICap + If + ICaT ) / Cm);
} // end CNZ




void NCZ_ODE_Initialise(double *state, int celltype) {
	state[0]  = -76.89967256;   // -77.13255836;
	state[1]  = 0.00584811327;   // 0.005631819916;
	state[2]  = 0.9126664364;   // 0.9168420281;
	state[3]  = 0.9339826669;   // 0.9380183564;
	state[4]  = 0.0002334162567;   // 0.0002267161277;
	state[5]  = 0.945924495;   // 0.9354212881;
	state[6]  = 0.001178151694;   // 0.00157418133;
	state[7]  = 0.02246300261;   // 0.02225979641;
	state[8]  = 9.838062008;   // 9.71;//10.30397012;
	state[9]  = 0.0001469561904;   // 0.0001403133065;
	state[10] = 140;   // 140;//140;//131.867138;
	state[11] = 0.7139970002;   // 0.7270823666;
	state[12] = 0.000839574626;   // 0.01223706011;
	state[13] = 0.9594992808;   // 0.8849139842;
	state[14] = 0.000262;   // 0.000262;
	state[15] = 0.948834;   // 0.948834;
	state[16] = 0.0001401836151;   // 0.0001313595105;
	state[17] = 0.9795264332;   // 0.9892411621;
	state[18] = 0.966862845;   // 0.9779168037;
	state[19] = 0.009453799877;   // 0.00958584702;
	state[20] = 0.009248109209;   // 0.009386941118;
	state[21] = 0.0005286378354;   // 0.000456694441;
	state[22] = 0.9350007301;   // 0.9571976502;
	state[23] = 0.1463722498;   // 0.1419271573;
	state[24] = 0.0003519550482;   // 0.0003667999273;
	state[25] = 0.9811908255;   // 0.9790238698;
	state[26] = 0.3122778037;   // 0.2977219443;
	state[27] = 0.002726913902;   // 0.002726913902;
	state[28] = 0.6643993051;   // 0.6643993051;
	state[29] = 0.0002484739247;   // 0.0002417881801;
	state[30] = 0.996802317;   // 0.9517278864;
	state[31] = 0.1992260523;   // 0.2000071464;
	state[32] = 0;   // 0;
	state[33] = 0;   // 0;
	state[34] = 1.901055114e-16;   // 2.355477517e-16;
	state[35] = 3.256315652e-15;   // 3.062722822e-15;
	state[36] = 0.9693633469;   // 0.9517278864;
}
