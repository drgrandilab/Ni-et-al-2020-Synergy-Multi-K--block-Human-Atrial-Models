/*

 This is the updatd CNZ model, update from the Colman et al. model 2013 Journal of Physiology Paper,

 Updates:  New IKur Model,
 Modifications to the conductances of several key ionic currents.

 Modifications to the Calcium dynamics model. change in the Refractory Parameters to make it faster, and thus
 remove the alternans seen at 400 ms of BCL.

 RMP: around -76 mV.

*/




#include "Updated_CNZ.h"
#include "ExplicitSolver.hpp"
#define CNZ_CON_STATE_LENGTH (34)



double Updated_CNZ_ODE(double dt, double Istim, SingleCellPara & het, double *state) {

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

	double Ena, Ek, Eca, Ekf, Enaf;
	double alpha, beta, tau, inf, a, b;
	double INa, IKr, IKs, ICaL, IK1, Iab, IbK, IbCa;
	double IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT;
	double fnak, sigma;
	double naidot, kidot, caidot;


	// Reversal Potentials
	Ena = 26.71 * log(CRN_nac / nai);
	Ek = 26.71 * log(CRN_kc / ki);
	Eca = 13.35 * log(CRN_cac / cai);

	//INa
	INa = het.GNa * Cm * CRN_gna * m * m * m * h * j * (V - Ena);

	// m gate
	alpha =  0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13)));
	if (fabs(V + 47.13) < 1e-10) alpha = 3.2;
	beta = 0.08 * exp(-V / 11.0);

	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	m = inf + (m - inf) * exp(-dt / 2.0 / tau); // steady state approx due to fast tau
	// m = inf + (m - inf) * exp(-dt / tau); // steady state approx due to fast tau
	m = inf + (m - inf) * exp(-dt / 2.0 / tau); // steady state approx due to fast tau
	//m = table_row[INA_M_INF] + (m - table_row[INA_M_INF])*table_row[INA_M_TAU];

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
	// IKr = 0.9*het.GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4));// 90% haibo
	IKr = 0.9 * het.GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4)); // 90% haibo

	// xr gate
	a = 0.0003 * (V + 14.1) / (1 - exp((V + 14.1) / -5));
	b = 0.000073898 * (V - 3.3328) / (exp((V - 3.3328) / 5.1237) - 1);
	if (fabs(V + 14.1) < 1e-10) a = 0.0015; /* denominator = 0 */
	if (fabs(V - 3.3328) < 1e-10) b = 3.7836118e-4; /* denominator = 0 */
	tau = 1 / (a + b);
	inf = 1 / (1 + exp((V + het.IKr_ac_shift + 14.1) / (-6.5 * het.IKr_grad)));
	xr = inf + (xr - inf) * exp(-dt / tau);

	// IKs
	IKs = het.GKs * Cm * CRN_gks * xs * xs * (V - Ek); // 1.1 times larger. haibo

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
	ICaL = 1.75 * het.GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL);   // haibo 1.75

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
	Ito = het.Gto * Cm * MAL_gto * itr * its * (V - Ek);

	//r gate
	inf = 1.0 / (1.0 + exp((V - 1.0) / -11.0));
	tau = (0.0035 * exp(-(V / 30.0) * 2) + 0.0015);
	itr = inf + (itr - inf) * exp(-(dt / 1000) / tau);
	//itr = table_row[ITO_ITR_INF] + (itr - table_row[ITO_ITR_INF])*table_row[ITO_ITR_TAU];

	// s gate
	inf = 1.0 / (1.0 + exp((V + 40.5) / 11.5));
	tau = (0.025635 * exp (-((V + 52.45) / 15.8827) * 2) + 0.01414);
	its = inf + (its - inf) * exp(-(dt / 1000) / tau);
	//its = table_row[ITO_ITS_INF] + (its - table_row[ITO_ITS_INF])*table_row[ITO_ITS_TAU];

	// 	//IKur
	// 	double IKur_K0 = -8.26597;
	// double IKur_c = 4.5128;
	// double IKur_x0 = 1.899769;
	// double IKur_y0 = 20.5232;

	IKur = het.GKur * Cm * CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597)))) * CNZ_a * CNZ_i * (V - Ek);
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


	inf = 1 / (1 + exp(-(V  + 5.52) / (8.81))); // change to simple version, haibo.  removed IKur mutations here.
	tau = (45.6666746826 / (1 + exp((V + 11.2306497073) / 11.5254705962)) + 4.26753514993) * (0.262186042981 / (1 + exp((V + 35.8658312707) / (-3.87510627762))) + 0.291755017928); //
	tau = tau / K_Q10;


	CNZ_a = inf + (CNZ_a - inf) * exp(-(dt) / tau);

	// CNZ_i
	inf = (het.IKur_inac_mult * 0.52424) / (1.0 + exp((V + 15.1142 + het.IKur_inac_shift) / (7.567021 * het.IKur_inac_grad))) + 0.4580778 + het.IKur_inac_add;
	tau = 2328 / (1 + exp(((V) - 9.435) / (3.5827))) + 1739.139;
	tau = tau / K_Q10;


	CNZ_i = inf + (CNZ_i - inf) * exp(-(dt) / tau);



	//If
	If = 0 * het.Gf * Cm * CRN_gf * If_y * (V - (-22));

	//If y gate
	inf = 1 / (1 + exp((V + 90.95 + het.If_vshift) / (10.1 * het.If_grad))   );
	tau = 1 / (1.2783E-04 * exp(-V / 9.2424) + 121.6092 * exp(V / 9.2424)   );
	If_y = inf + (If_y - inf) * exp(-dt / tau);

	//ICaT
	ICaT = het.GCaT * Cm * gcaT * ff * dd * (V - EcaT);

	// dd gate
	a = 1.0680 * exp((V + 26.3) / 30.0);
	b  = 1.0680 * exp(-(V + 26.3) / 30.0);
	tau = 1.0 / (a + b);
	inf = 1.0 / (1.0 + exp(-(V + 37.0) / 6.8));
	dd = inf + (dd - inf) * exp(-dt / tau);

	// ff gate
	a = 0.0153 * exp(-(V + 71.0) / 83.3);
	b  = 0.0150 * exp((V + 71.0) / 15.38);
	tau = 1.0 / (a + b);
	inf = 1.0 / (1.0 + exp((V + 71.0) / 9.0));
	ff = inf + (ff - inf) * exp(-dt / tau);

	//IK1
	IK1 = het.GK1 * Cm * CRN_gk1 * (V - Ek + het.IK1_v_shift) / (1 + exp(0.07 * (V + 80 + het.IK1_v_shift)));
	//IK1 = GK1 * Cm * CRN_gk1 * (V-Ek+het.IK1_v_shift) * (het.IK1_type_default * (table_row[IK1_FAC]) + het.IK1_type_AS *(table_row[IK1_FAC_AS])  );

	//Iab
	Iab = Cm * 0.0003879 * (V + 69.6) / (1 - 0.8377 * exp((V + 49.06) / 1056));

	//Ibk
	IbK =  0.0;//Cm * CRN_gbk * (V - Ek);

	//IbCa
	IbCa = 1.05 * Cm * CRN_gbca * (V - Eca); // increased background IbCa, haibo
	// IbNa
	IbNa = 1.5 * Cm * CRN_gbna * (V - Ena);  // increase. haibo

	//ICap
	ICap = 1.0 * het.GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);  // change, haibo
	// ICap = 1.1 * het.GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);  // change, haibo


	//INaCa, times 1.5, Haibo
	INaCa = 1.5 * het.GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
	        (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
	        (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
	         exp(V * (CRN_gammalr - 1) * F / (R * T)));


	// Version 3.0
	/*INaCa = 1.6 * het.GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
	        (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
	        (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
	         exp(V * (CRN_gammalr - 1) * F / (R * T)));*/


	//INaK
	sigma = (exp(CRN_nac / 67.3) - 1) / 7.0;
	fnak = 1 / (1 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0365 * sigma * exp(-V * F / (R * T)));
	// INaK = 1.3 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(/*CRN_kmnai*/11/ nai, 4));  // increase, Haibo.

	// verison 3.0
	INaK = 1.3 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(/*CRN_kmnai*/10.5 / nai, 4)); // increase, Haibo.

	// INaK = 1.28 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(CRN_kmnai / nai, 1.5));

	// concentrations
	naidot = (-3 * INaK - 3 * INaCa - IbNa - INa) / (F * CRN_vi);
	kidot = (2 * INaK - IK1 - Ito - IKur - IKr - IKs - IbK) / (F * CRN_vi);
	nai = nai + dt * naidot;
	// ki = ki + dt * kidot;

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
	Jj_nj = 2.1 * DCa * Aj_nj / xj_nj * ((Cass) - cai) * 1e-6;  // increase, Haibo

	double J_SERCASR, J_bulkSERCA;
	double J_SERCASRss, J_bulkSERCAss;

	J_SERCASR =  0.7 * (-k3 * pow(CaSR1, 2) * (cpumps - SERCACa) + k4 * (SERCACa)) * Vnonjunct3 * 2;
	J_bulkSERCA =   0.7 * (k1 * pow(cai, 2) * (cpumps - (SERCACa)) - k2 * (SERCACa)) * Vnonjunct3 * 2;   // increase, haibo

	J_SERCASRss = 0.7 * (-k3 * pow(CaSR2, 2) * (cpumps - (SERCACass)) + k4 * (SERCACass)) * Vss * 2;
	J_bulkSERCAss = 0.7 * (k1 * pow((Cass), 2) * (cpumps - (SERCACass)) - k2 * (SERCACass)) * Vss * 2;   // increase, haibo



	// change parameters, to remove alternans.
	double RyRtauadapt = 0.3;
	double RyRtauactss = 5e-3;
	double RyRtauinactss = 1 * 15e-3;
	double RyRtauact =  5e-3;//0.6 * 8.75e-3;
	double RyRtauinact = 2 * 15e-3; //0.6 * 87.5e-3;

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
	state[33] = If_y;
	return (Istim + (INa + Ito + IKur + IKr + IKs + ICaL + IK1 + IbK + IbNa + IbCa + Iab + INaCa + INaK + ICap + If + ICaT ) / Cm);
} // end CNZ


double Updated_CNZ_ODE_Temp(double dt, double Istim, SingleCellPara & het, double *state) {

	double Ena, Ek, Eca, Ekf, Enaf;
	double alpha, beta, tau, inf, a, b;
	double INa, IKr, IKs, ICaL, IK1, Iab, IbK, IbCa;
	double IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT;
	double fnak, sigma;
	double naidot, kidot, caidot;

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
	double If_y = state[33];
	// double rkv  = state[29];
	// double skv  = state[30];
	// double kif  = state[31];
	// double naif = state[32];
	// double Vmf  = state[33];
	// double CNZ_a = state[34];
	// double CNZ_i = state[35];

	// double If_y = state[38];

	// Reversal Potentials
	Ena = 26.71 * log(CRN_nac / nai);
	Ek = 26.71 * log(CRN_kc / ki);
	Eca = 13.35 * log(CRN_cac / cai);
	// Ekf = 26.71 * log(kof / kif);
	// Enaf = 26.71 * log(naof / naif);

	//INa
	INa = het.GNa * Cm * CRN_gna * m * m * m * h * j * (V - Ena);

	// m gate
	alpha =  0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13)));
	if (fabs(V + 47.13) < 1e-10) alpha = 3.2;
	beta = 0.08 * exp(-V / 11.0);

	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	// m = inf + (m - inf) * exp(-dt / tau); // steady state approx due to fast tau
	m = inf + (m - inf) * exp(-dt / 2.0 / tau); // steady state approx due to fast tau   // to make it more accurate ...
	m = inf + (m - inf) * exp(-dt / 2.0 / tau); // steady state approx due to fast tau

	// j gate
	if (V >= -40.0) {
		alpha  = 0.0;
		beta = 0.3 * exp(-2.535e-7 * V) / (1 + exp(-0.1 * (V + 32)));
	} else {
		alpha = (-1.2714e5 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) * (V + 37.78) / (1 + exp(0.311 * (V + 79.23)));
		beta = 0.1212 * exp(-0.01052 * V) / (1 + exp(-0.1378 * (V + 40.14)));
	}
	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	j = inf + (j - inf) * exp(-dt / tau);

	// h gate
	if (V >= -40.0) {
		alpha  = 0.0;
		beta = 1 / (0.13 * (1 + exp((V + 10.66) / -11.1)));
	} else {
		alpha = 0.135 * exp((V + 80) / -6.8);
		beta = 3.56 * exp(0.079 * V) + 3.1e5 * exp(0.35 * V);
	}
	tau = 1 / (alpha + beta);
	inf = alpha * tau;

	h = inf + (h - inf) * exp(-dt / tau);

	// IKr



	// wettwer . et al. 2004
	/* double fvm = 1.0 / (1.0 + exp((V + 9.0) / 22.4));
	 double Wettwer_gkr = 0.05; // ms/uF
	 Wettwer_gkr = CRN_gkr; // ms/uF
	 // xr = state[6];
	 // xr gate
	 double vshift = 0;
	 if (fabs(V - vshift + 14.2) < 1e-10) a = 0.0015; // denominator = 0
	 else
	     a = 0.00138 * (V - vshift + 14.2) / (1.0 - exp((V - vshift + 14.2) / -8.13));
	 if (fabs(V - vshift + 38.9) < 1e-10) b = 3.7836118e-4;  //denominator = 0
	 else
	     b = 0.00061 * (V - vshift + 38.9) / (exp((V - vshift + 38.9) / 6.9) - 1.0);

	 tau = 1 / (a + b);
	 inf = 1 / (1 + exp((V + IKr_ac_shift + 21.5) / (-7.5 * IKr_grad)));
	 // xr = inf + (xr-inf)*exp(-dt/tau);
	 // state[6] = xr;

	 xr = inf + (xr - inf) * exp(-dt / tau);
	 // IKr = gkr * Y[14] * rkr * (V - EK);

	 IKr = 0.5*GKr * Cm * Wettwer_gkr * fvm * xr * (V - Ek);*/





	// Grandi IKr
	/*    double GB_Ko = CRN_kc;
	    double gkr = 0.035 * sqrt(GB_Ko / 5.4);
	    inf = 1.0 / (1.0 + exp(-(V + 10.0) / 5.0));
	    tau = 550.0 / (1.0 + exp((-22.0 - V) / 9.0)) * 6.0 / (1.0 + exp((V - (-11.0)) / 9.0))
	          + 230.0 / (1.0 + exp((V - (-40.0)) / 20.0));


	    xr = inf + (xr - inf) * exp(-dt / tau);
	    double rkr = 1.0 / (1.0 + exp((V + 74.0) / 24.0));

	    IKr = GKr * Cm * gkr * xr * rkr * (V - Ek);
	*/



	// CRN IKr

	// IKr = GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4));
	IKr = 0.9 * het.GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 9) / 22.4));
	a = 0.0003 * (V + 14.1) / (1 - exp((V + 14.1) / -5));
	b = 0.000073898 * (V - 3.3328) / (exp((V - 3.3328) / 5.1237) - 1);
	if (fabs(V + 14.1) < 1e-10) a = 0.0015;
	if (fabs(V - 3.3328) < 1e-10) b = 3.7836118e-4;  //
	tau = 1 / (a + b);
	inf = 1 / (1 + exp((V + het.IKr_ac_shift + 14.1) / (-6.5 * het.IKr_grad)));
	xr = inf + (xr - inf) * exp(-dt / tau);
	// end CRN IKr
	// IKs

	// Ek = 26.71 * log(CRN_kc / ki);

	double pNaK = 0.01833;
	double eks = 26.71 * log((CRN_kc + pNaK * CRN_nac) / (ki + pNaK * nai));
	// IKs = GKs * Cm * CRN_gks * xs * xs * (V - eks);  // use grandi style of reversal potential
	IKs = 1.1 * het.GKs * Cm * CRN_gks * xs * xs * (V - Ek);

	// xs gate
	a = 0.00004 * (V - 19.9) / (1 - exp((V - 19.9) / -17));
	b = 0.000035 * (V - 19.9) / (exp((V - 19.9) / 9) - 1);
	if (fabs(V - 19.9) < 1e-10) { /* denominator = 0 */
		a = 0.00068;
		b = 0.000315;
	}
	tau = 0.5 / (a + b); // note lagrer taus may be more accurate
	inf = sqrt(1 / (1 + exp((V - 19.9 - het.IKs_shift) / (-12.7 * het.IKs_grad))));
	xs = inf + (xs - inf) * exp(-dt / tau);

	//ICaL
	ICaL = /*2.125 **/ 1.75 * het.GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL);

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
	Ito = het.Gto * Cm * MAL_gto * itr * its * (V - Ek);

	// r gate
	inf = 1.0 / (1.0 + exp((V - 1.0) / -11.0));
	tau = (0.0035 * exp(-(V / 30.0) * 2) + 0.0015);
	itr = inf + (itr - inf) * exp(-(dt / 1000) / tau);

	// s gate
	inf = 1.0 / (1.0 + exp((V + 40.5) / 11.5));
	tau = (0.025635 * exp (-((V + 52.45) / 15.8827) * 2) + 0.01414);
	its = inf + (its - inf) * exp(-(dt / 1000) / tau);

	//IKur
	// IKur_type_CNZ = 0;

	double IKur_K0 = -8.26597;
	double IKur_c = 4.5128;
	double IKur_x0 = 1.899769;
	double IKur_y0 = 20.5232;
	int IKur_type_CNZ = 1;
	IKur = ((1 - IKur_type_CNZ) * ( Cm * MAL_gkur * isusr * isuss * (V - Ek))) +
	       (IKur_type_CNZ) * Cm * CNZ_gkur
	       * /*GKur_scale * */
	       (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597))))
	       * CNZ_a * CNZ_i * (V - Ek);
	IKur = het.GKur * IKur;
	// printf("GKur = %f\n", GKur);
	// printf("CNZ_conductance_scaleF = %f\n", CNZ_conductance_scaleF);

	//isusr
	/*inf = 1 / (1 + exp(-(V + IKur_ac_shift + 6) / (8.6 * IKur_ac_grad)));
	tau = (0.009 / (1.0 + exp((V + 5) / 12)) + 0.0005);
	isusr = inf + (isusr - inf) * exp(-(dt / 1000) / tau);

	//isuss
	inf = 1 / (1 + exp((V + IKur_inac_shift + 7.5) / (10 * IKur_inac_grad)));
	tau = (0.59 / (1 + exp((V + 60.0) / 10)) + 3.05);
	isuss = inf + (isuss - inf) * exp(-(dt / 1000) / tau);*/


	double K_Q10 = 3.0;// 3.5308257834747638;
	//CNZ_a

	// inf = ((IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + IKur_ac1_shift) / (-5.75418 * IKur_ac1_grad))) ) * ((IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + IKur_ac2_shift) / (-11.51037561 * IKur_ac2_grad)))) + IKur_ac_add;
	// inf = 1 / (1 + exp(-(V + IKur_ac_shift + 5.52) / (8.81 * IKur_ac_grad)));
	inf = 1 / (1 + exp(-(V  + 5.52) / (8.81 )));


	tau = (45.6666746826 / (1 + exp((V + 11.2306497073) / 11.5254705962)) + 4.26753514993) * (0.262186042981 / (1 + exp((V + 35.8658312707) / (-3.87510627762))) + 0.291755017928); //
	tau = tau / K_Q10;

	CNZ_a = inf + (CNZ_a - inf) * exp(-(dt) / tau);

	//CNZ_i
	inf = (het.IKur_inac_mult * 0.52424) / (1.0 + exp((V + 15.1142 + het.IKur_inac_shift) / (7.567021 * het.IKur_inac_grad))) + 0.4580778 + het.IKur_inac_add;
	// inf = Tau_scale * (IKur_inac_mult * 0.52424) / (1.0 + exp((V + 15.1142 ) / (7.567021))) + 0.4580778 + IKur_inac_add;
	tau = 2328 / (1 + exp(((V) - 9.435) / (3.5827))) + 1739.139;
	/*if (V<-30)
	{
	    tau = 500.0;
	}*/


	tau = tau / K_Q10;

	CNZ_i = inf + (CNZ_i - inf) * exp(-(dt) / tau);

	/*// If  // No If here.
	If = 0.0 * het.Gf * Cm * CRN_gf * If_y * (V - (-22));

	//If y gate
	inf = 1 / (1 + exp((V + 90.95 + het.If_vshift) / (10.1 * het.If_grad))   );
	tau = 1 / (1.2783E-04 * exp(-V / 9.2424) + 121.6092 * exp(V / 9.2424)   );

	If_y = inf + (If_y - inf) * exp(-dt / tau);*/


	If = 0.0 ;//* het.Gf * Cm * CRN_gf * If_y * (V - (-22));

	//ICaT
	ICaT = het.GCaT * Cm * gcaT * ff * dd * (V - EcaT);

	// dd gate
	a = 1.0680 * exp((V + 26.3) / 30.0);
	b  = 1.0680 * exp(-(V + 26.3) / 30.0);
	tau = 1.0 / (a + b);
	inf = 1.0 / (1.0 + exp(-(V + 37.0) / 6.8));
	dd = inf + (dd - inf) * exp(-dt / tau);

	// ff gate
	a = 0.0153 * exp(-(V + 71.0) / 83.3);
	b  = 0.0150 * exp((V + 71.0) / 15.38);
	tau = 1.0 / (a + b);
	inf = 1.0 / (1.0 + exp((V + 71.0) / 9.0));
	ff = inf + (ff - inf) * exp(-dt / tau);

	//IK1
	IK1 = het.GK1 * Cm * CRN_gk1 * (V - Ek + het.IK1_v_shift) / (1 + exp(0.07 * (V + 80 + het.IK1_v_shift)));

	//Iab
	Iab = 0.0 * Cm * 0.0003879 * (V + 69.6) / (1 - 0.8377 * exp((V + 49.06) / 1056));

	//Ibk
	IbK = 0.0 * Cm * CRN_gbk * (V - Ek);

	//IbCa
	IbCa = 1.05 * het.GbCa * Cm * CRN_gbca * (V - Eca);

	// IbNa
	IbNa = 1.5 * Cm * CRN_gbna * (V - Ena);

	//ICap
	ICap = 1.2 * het.GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);

	//INaCa
	INaCa = 1.5 * het.GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
	        (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
	        (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
	         exp(V * (CRN_gammalr - 1) * F / (R * T)));

	//INaK
	sigma = (exp(CRN_nac / 67.3) - 1) / 7.0;
	fnak = 1 / (1 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0365 * sigma * exp(-V * F / (R * T)));
	INaK = 1.28 * Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(CRN_kmnai / nai, 1.5));

	double IGap;
	IGap = 0;

	// Concentrations
	naidot = (-3 * INaK - 3 * INaCa - IbNa - INa) / (F * CRN_vi);
	kidot = (2 * INaK - IK1 - Ito - IKur - IKr - IKs - IbK) / (F * CRN_vi);
	nai = nai + dt * naidot;
	ki = ki + dt * kidot;

	//Caclium handling
	double betass;
	betass = pow(( 1 + SLlow * KdSLlow / pow(((Cass) + KdSLlow), 2) + SLhigh * KdSLhigh / pow(((Cass) + KdSLhigh), 2) + BCa * KdBCa / pow(((Cass) + KdBCa), 2)  ), (-1));

	double betai, gammai;
	betai = pow(( 1 + BCa * KdBCa / pow((cai + KdBCa), 2)  ), (-1));
	gammai = BCa * KdBCa / pow((cai + KdBCa), 2);

	double betaSR1, betaSR2;
	betaSR1 = pow( ( 1 + CSQN * KdCSQN / pow((CaSR1 + KdCSQN), 2) ), (-1));
	betaSR2 = pow( ( 1 + CSQN * KdCSQN / pow((CaSR2 + KdCSQN), 2) ), (-1));

	double Jj_nj;
	Jj_nj = 2.1 * DCa * Aj_nj / xj_nj * ((Cass) - cai) * 1e-6;

	double J_SERCASR, J_bulkSERCA;
	double J_SERCASRss, J_bulkSERCAss;

	J_SERCASR =  (-k3 * pow(CaSR1, 2) * (cpumps - SERCACa) + k4 * (SERCACa)) * Vnonjunct3 * 2;
	J_bulkSERCA =   1.2 * (k1 * pow(cai, 2) * (cpumps - (SERCACa)) - k2 * (SERCACa)) * Vnonjunct3 * 2;

	J_SERCASRss = (-k3 * pow(CaSR2, 2) * (cpumps - (SERCACass)) + k4 * (SERCACass)) * Vss * 2;
	J_bulkSERCAss = 1.2 * (k1 * pow((Cass), 2) * (cpumps - (SERCACass)) - k2 * (SERCACass)) * Vss * 2;

	/*double Q10SRCaP = 2.6;
	double Vmax_SRCaP = 5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
	double Kmf = 2.5*0.246e-3;          // [mM] was    Kmf = 2.0*0.246e-3; // Haibo
	double Kmr = 1.7;               // [mM]L cytosol
	double hillSRCaP = 1.787;       // [mM]
	J_bulkSERCA= 0.5*pow(Q10SRCaP, 1.0) * Vmax_SRCaP * (pow(cai/ Kmf, hillSRCaP) - pow(CaSR1/ Kmr, hillSRCaP)) / (1.0 + pow(cai / Kmf, hillSRCaP) + pow(CaSR1 / Kmr, hillSRCaP));
	J_bulkSERCAss=0.5*pow(Q10SRCaP, 1.0) * Vmax_SRCaP * (pow(Cass/ Kmf, hillSRCaP) - pow(CaSR2/ Kmr, hillSRCaP)) / (1.0 + pow(Cass / Kmf, hillSRCaP) + pow(CaSR2 / Kmr, hillSRCaP));*/

	double RyRtauadapt = 0.25 * 1;
	double RyRtauactss = 5e-3;
	double RyRtauinactss = 15e-3;
	double RyRtauact = 0.7 * 8.75e-3;
	double RyRtauinact = 0.7 * 87.5e-3;


	// original:
	/* double RyRtauadapt = 1;

	 double RyRtauactss = 5e-3;
	 double RyRtauinactss = 15e-3;
	 double  RyRtauact = 8.75e-3;
	 double RyRtauinact = 87.5e-3;*/

	/*
	    double RyRtauadapt = 0.25*1;

	    double RyRtauactss = 0.5*5e-3;
	    double RyRtauinactss = 0.5*15e-3;
	    double  RyRtauact = 0.8*18.75e-3;
	    double RyRtauinact =0.8* 87.5e-3;*/

	double nuss = 1 * 625 * Vss;
	double RyRSRCass = (1 - 1 / (1 +  exp((CaSR2 - 0.3) / 0.1)));
	double RyRainfss = 0.505 - 0.427 / (1 + exp((( Cass + (het.fRyR * Cass) ) * 1000 - 0.29) / 0.082));
	double RyRoinfss = (1 - 1 / (1 +  exp(((Cass + (het.fRyR * Cass)   ) * 1000 - ((RyRass) + 0.22)) / 0.03)));
	double RyRcinfss = (1 / (1 + exp(((Cass + (het.fRyR * Cass )) * 1000 - ((RyRass) + 0.02)) / 0.01)));
	double Jrelss = nuss * ( (RyRoss) ) * (RyRcss) * RyRSRCass * ( CaSR2 -  (Cass) );

	double nu3 = 1 * Vnonjunct3;
	double RyRSRCa3 = (1 - 1 / (1 +  exp((CaSR1 - 0.30) / 0.1)));
	double RyRainf3 =  0.505 - 0.427 / (1 + exp(( (cai * 2 + ( het.fRyR * cai * 2)  ) * 1000 - 0.29) / 0.082));
	double RyRoinf3 = (1 - 1 / (1 +  exp(( (cai * 2 + ( het.fRyR * cai * 2) ) * 1000 - ((RyRa3) + 0.22)) / 0.03)));
	double RyRcinf3 = (1 / (1 +  exp(( (cai * 2 + (het.fRyR * cai * 2 ) ) * 1000 - ((RyRa3) + 0.02)) / 0.01)));
	double Jrel3 = nu3 * ( (RyRo3) ) * (RyRc3) * RyRSRCa3 * ( CaSR1 -  cai );

	Jrelss = het.fIRel * Jrelss;
	Jrel3 = het.fIRel * Jrel3;

	double JSRCaleak3 = 0.5 * het.GSR_leak * kSRleak * ( CaSR1 - cai ) * Vnonjunct3;
	double JSRCaleakss = 0.5 * het.GSR_leak * kSRleak * ( CaSR2 - (Cass) ) * Vss;

	// O_JSERCASR = J_SERCASR;
	// O_JSERCASRss = J_SERCASRss;
	// O_JBULKSERCA = J_bulkSERCA;
	// O_JBULKSERCAss = J_bulkSERCAss;
	// O_JRELss = Jrelss;
	// O_JREL = Jrel3;
	// O_JSRLEAK = JSRCaleak3;
	// O_JSRLEAKss = JSRCaleakss;

	double JCa, JCass;
	JCa = -/*BULK_CONST **/ J_bulkSERCA + JSRCaleak3 + Jrel3 + Jj_nj;
	JCass = -Jj_nj + JSRCaleakss - /*BULK_CONST * */J_bulkSERCAss + Jrelss;

	double JSRCa1, JSRCa2;
	JSRCa1 = J_SERCASR - JSRCaleak3 - Jrel3;
	JSRCa2 = J_SERCASRss - JSRCaleakss - Jrelss;

	double dy;

	dy = /*0.5 * */(-J_SERCASR + J_bulkSERCA) / Vnonjunct3;
	SERCACa = Foward_Euler(dt / 1000, dy, SERCACa);
	dy = /*0.5 * */(-J_SERCASRss + J_bulkSERCAss) / Vss;
	SERCACass = Foward_Euler(dt / 1000, dy, SERCACass);

	RyRoss = Euler_inf(dt / 1000, RyRoss, RyRoinfss, RyRtauactss);
	RyRcss = Euler_inf(dt / 1000, RyRcss, RyRcinfss, RyRtauinactss);
	RyRass = Euler_inf(dt / 1000, RyRass, RyRainfss, RyRtauadapt);
	RyRo3 = Euler_inf(dt / 1000, RyRo3, RyRoinf3, RyRtauact);
	RyRc3 = Euler_inf(dt / 1000, RyRc3, RyRcinf3, RyRtauinact);
	RyRa3 = Euler_inf(dt / 1000, RyRa3, RyRainf3, RyRtauadapt);

	dy =  betass * ( JCass / Vss + ((-( /*RYR **/ ICaL) - IbCa - ICap - ICaT + 2 * INaCa)) / (2 * Vss * 1000 * F) );
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

	// state[29] = rkv;
	// state[30] = skv;
	// state[31] = kif;
	// state[32] = naif;
	// state[33] = Vmf;
	state[29] = CNZ_a;
	state[30] = CNZ_i;

	state[33] = If_y;
// (Istim + (INa + Ito + IKur + IKr + IKs + ICaL + IK1 + IbK + IbNa + IbCa + Iab + INaCa + INaK + ICap + If + ICaT ) / Cm);
	double Iion_tot = Istim + (INa + Ito + IKur + IKr + IKs + ICaL + IK1 + IbK + IbNa + IbCa + Iab + INaCa + INaK + ICap + If + ICaT /*+ (FB_number * IGap)*/) / Cm;
	// (INa + Ito + IKur + IKr + IKs + ICaL + IK1 + IbK + IbNa + IbCa + Iab + INaCa + INaK + ICap)/Cm;
	// current_out = Ito / Cm;
	return Iion_tot;

} // end Return Iion_tot




void Updated_CNZ_ODE_Initialise(double *state, int celltype) {
	// state[0] = -76.079842; // V
	// state[1] = 0.006676; // m
	// state[2] = 0.896736; // h
	// state[3] = 0.918836; // j
	// state[4] = 0.000259; // d
	// state[5] = 0.999059; // f
	// state[6] = 0.000072; // xr
	// state[7] = 0.022846; // xs
	// state[8] = 11.170000; // nai
	// state[9] = 0.000098;  // Cai
	// state[10] = 115.632438; // ki
	// state[11] = 0.770253; // fca
	// state[12] = 0.000905; // Itr
	// state[13] = 0.956638; // Its
	// state[14] = 0.000289; // Isusr
	// state[15] = 0.998950; // Isuss
	// state[16] = 0.000104; // Cass
	// state[17] = 0.437859; // CaSR1
	// state[18] = 0.434244; // CaSR2
	// state[19] = 0.002432; // SERCACa
	// state[20] = 0.002443; // SERCACass
	// state[21] = 0.000412; // RyRoss
	// state[22] = 0.967156; // RyRcss
	// state[23] = 0.118222; // RyRass
	// state[24] = 0.000365; // RyRo3
	// state[25] = 0.977008; // RyRc3
	// state[26] = 0.115451; // RyRa3
	// state[27] = 0.003182; // dd
	// state[28] = 0.637475; // ff
	// state[29] = 0; // CNZ_a
	// state[30] = 1;//0.9; // CNZ_i
	// state[31] = 0; // if y
	// state[32] = 0; //ACh_i
	// state[33] = 0; // ACh_j

	state[0]  = -76.43581147;
	state[1]  = 0.006303231435;
	state[2]  = 0.9040359445;
	state[3]  = 0.9261072408;
	state[4]  = 0.0002473410325;
	state[5]  = 0.9369631854;
	state[6]  = 0.001838355913;
	state[7]  = 0.02287634678;
	state[8]  = 11.28749024;
	state[9]  = 0.0001648565885;
	state[10] = 131.867138;
	state[11] = 0.6975282552;
	state[12] = 0.01160697603;
	state[13] = 0.8955390741;
	state[14] = 0.000262;
	state[15] = 0.948834;
	state[16] = 0.0001517585657;
	state[17] = 1.185296236;
	state[18] = 1.171510414;
	state[19] = 0.0124452557;
	state[20] = 0.01220050338;
	state[21] = 0.0004635804651;
	state[22] = 0.9560508579;
	state[23] = 0.1617850925;
	state[24] = 0.0001972118928;
	state[25] = 0.9966475449;
	state[26] = 0.3654457351;
	state[27] = 0.003020243476;
	state[28] = 0.6469698391;
	state[29] = 0.00031910821;
	state[30] = 0.9499663408;
	state[31] = 0.199622;
	state[32] = 0.000000;
	state[33] = 0.000000;



}


int cnz_con_length() {
	return CNZ_CON_STATE_LENGTH;
}

