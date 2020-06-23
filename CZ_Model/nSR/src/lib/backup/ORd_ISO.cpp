#include "ORd_ISO.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define  zca (2.0)
#define zna (1.0)

double INa_Drug_ORd_Model(double dt, double stims, double *state, TypeCell celltype, int ISO, IKsMutation mutation, float ABIndex,  SingleCellPara& het)
{
	double ENa, EK, EKs;
	double INa, INaL, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCa_i, INaCa_ss, INaCa, INaK, IKb, INab, IpCa, ICab, Ist;
	double Jrel, Jup, Jtr, Jdiff, JdiffNa, JdiffK, Jleak;
	double CaMKa, CaMKb;
	Ist = stims;

	double INah_half, INa_CaMh_half, GNa_ISO;

	double ICaLd_half, ICaLf_half, Pca_ISO;
	double Tau_ICaL_ff_ISO;

	double Tau_Ks, GKs_ISO, IKs_ac_shift_ISO;

	double NaK_ISO;

	double GKb_ISO;


	double A_rel_ISO, tau_rel_ISO;
	double Jup_ISO;
	double kmtrpn_ISO;

	INah_half = INa_CaMh_half = 0.0;
	GNa_ISO = 1.0;
	ICaLd_half = 0.0;
	ICaLf_half = 0.0;
	Tau_ICaL_ff_ISO = 1.0;
	Pca_ISO = 1.0;
	GKs_ISO = 1.0;
	Tau_Ks = 1.0;
	NaK_ISO = 1.0;//
	GKb_ISO = 1.0;
	A_rel_ISO = 1.0;   // ref 2
	tau_rel_ISO = 1.0;  // ref 2
	Jup_ISO = 1.0;
	kmtrpn_ISO = 1.0;

	CaMKb = CaMKo * (1.0 - state[40]) / (1.0 + KmCaM / state[6]);
	CaMKa = CaMKb + state[40];
	state[40] += dt * (aCaMK * CaMKb * (CaMKb + state[40]) - bCaMK * state[40]);




	if (ISO == 1)
	{
		double time_1;
		/*        if (t > 100000)
		        {
		            time_1 = t - 100000;

		        }
		        else {
		            time_1 = 0.0;
		        }
		*/
		double value_1 = 1.0;// - exp(-time_1 / 7700);
		double value_2 = 1.0;// - exp(-time_1 / 39700);

		/*if (time_1 > 100000)
		{
		    value_2 = value_1 = 1.0;
		}*/
		INah_half = INa_CaMh_half = value_1 * 5.0;
		GNa_ISO =   1.0 + value_1 * 1.70;
		ICaLd_half = value_1 * 9.9; // ref 1 and 2 16.0
		// change to 8 will give good results  // the original value from rudy paper was 16.0
		/*Figure 3 in paper         // from  paper : Chen et al. L-Type Ca2+ Channel Density and Regulation Are Altered in Failing Human Ventricular Myocytes and Recover After Support With Mechanical Assist Devices
		ISO caused a significant leftward shift of the voltage dependence of I Ca,L activation in
		NF- (change of V 0.5 (dρ): -9.9 mV; P<0.05) and LVAD- (change of V 0.5 (dρ): -5.8 mV; P0.05) HVMs but no significant shift in F- (change of V 0.5 (dρ): -2.2 mV; PϾ0.05) HVMs.*/
		// td    = 0.6 + 1.0 / (exp(-0.05 * (v + 6.0 + 16.0)) + exp(0.09 * (v + 14.0 + 16.0)));  // in the rudy paper, the td was not changed by pka..
		Tau_ICaL_ff_ISO = 1 - value_1 * 0.5; // ISO accerartes the fast inactivation from paper chen et al. L-Type Ca2+ Channel Density and Regulation Are Altered in Failing Human Ventricular Myocytes and Recover After Support With Mechanical Assis
		/* In the paper */
		ICaLf_half = value_1 * 8.0; //ref 1
		Pca_ISO =  1.0 + value_1 * 1.5;

		GKs_ISO = 1.0 + value_2 * 2.2;
		IKs_ac_shift_ISO = value_2 * 7.0; // Impaired IKs channel activation by Ca2 +-dependent PKC shows correlation with emotion/arousal-triggered events in LQT1

		if (mutation == ORd or mutation == IKS_WT)
		{
			GKs_ISO = 1.0 + value_2 * 2.2;
			Tau_Ks = 1.0 - value_2 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
			IKs_ac_shift_ISO = 7.0 * value_2;
		} else if (mutation == HET)
		{
			GKs_ISO = 1 + value_2 * 0.88;
			Tau_Ks = 1.0 - value_2 * 0.4 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
			IKs_ac_shift_ISO = 7.0 * 0.4 * value_2;
			// GKs_ISO = 3.2 * 0.4;
		} else {
			// GKs_ISO = 3.2 * 0.4;
			// GKs_ISO = 2.23;
			// GKs_ISO = 1.0;
			GKs_ISO = 1 + value_2 * 0.88;
			Tau_Ks = 1.0 - value_2 * 0.4 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
			IKs_ac_shift_ISO = 7.0 * 0.4 * value_2;

		}
		// / say: previous work in human [(20) using original current formulations (21) with τxs1,PKA = 0.6 * τxs1 and GKs,PKA = 3.2 * GKs,
		NaK_ISO = 1 - value_1 * 0.3; //
		GKb_ISO = 1 + value_1 * 2.62; // ref 2 3.62 and ref 1 2.5
		A_rel_ISO = 1 + value_1 * 0.8; // ref 2
		tau_rel_ISO = 1 - value_1 * 0.6; // ref 2
		Jup_ISO = 1 - value_1 * 0.46;
		kmtrpn_ISO = 1.0 + value_1 * 0.6;
		/*GKs_ISO = 1.0 + value_2 * 2.2;
		Tau_Ks = 1.0 - value_2 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
		IKs_ac_shift_ISO = 7.0 * value_2;*/
	}


	/*  default values for ABIndex = 0.0 (Value set 0 == Apical Cells)
	        and ABIndex ranges from 0.0 -> 1.0 */
	/*
	*   0.0 -- Apical Cells
	*   1.0 -- Basal Cells
	*/

	double GTo_ABh = 1.0;
	double InAcTo_Vhalf_ABh = 0.0;
	double InAcTo_Vk_ABh = 1.0;
	double GKs_ABh = 1.0;
	double TauKs_ABh = 1.0;
	double AcKs_Vhalf_ABh = 0.0;
	double AB_IKr_ActVhalf = 0.0;




	if (ABIndex < 0.0 or ABIndex > 1.0) {
		std::cerr << "Wrong Apcial-Basal Ratio, at ABIndex = " << ABIndex << std:: endl;
		std::cerr << "Program Existing... " << std::endl;
		std::exit(0);
	}
	/* V1, ABIndex   == 0.0 -> Apical Cells */
	GTo_ABh          =  1.0 - ABIndex * (1 - (16.5) / 29.6) ; //(1.0 + (1 - ABIndex) * (/ 16.529.6  - 1));      // Apical cells: 1.0; Basal Cells: 16.5/29.6
	InAcTo_Vhalf_ABh = ABIndex * 4.0;                          // Apical Cells: 0; Basal Cells: 4.0;
	// InAcTo_Vk_ABh    = 1.0 + (1 - ABIndex) * ((4.5 - 3.4) / 3.4);   // Apical cells: 1.0; Basal Cells: 1.0  // not used here..
	GKs_ABh          = 1 - ABIndex * (1.0 - 2.1 / 5.6);  // Apical Cells: 1.0; Basal Cells: 2.1/5.6
	// TauKs_ABh        = 1.0 + (1 - ABIndex) * ((358.0 - 516.0) / 516.0);      //  // not used here.. refer to TNNP model

	// Iion_out *ion;
	// printf("%d\n",celltype);

	// double dt =  dti;
	// revpots(state);
	ENa = (R * T / F) * log(nao / state[1]);
	EK = (R * T / F) * log(ko / state[3]);
	EKs = (R * T / F) * log((ko + 0.01833 * nao) / (state[3] + 0.01833 * state[1]));

	// f// printf(stderr, "Ena=%f state10=%f\n", ENa, state[10]);
	// RGC(dt, state, celltype);
	CaMKb = CaMKo * (1.0 - state[40]) / (1.0 + KmCaM / state[6]);
	CaMKa = CaMKb + state[40];
	double vffrt = state[0] * F * F / (R * T);
	double vfrt = state[0] * F / (R * T);


	/* ORd INa*/

	/*double tm = 1.0 / (6.765 * exp((state[0] + 11.64) / 34.77) + 8.552 * exp(-(state[0] + 77.42) / 5.955));
	{   double mss = 1.0 / (1.0 + exp((-(state[0] + 39.57)) / 9.871));
	    state[9] = mss - (mss - state[9]) * exp(-dt / tm);
	    double hss = 1.0 / (1 + exp((state[0] + 82.90 + INah_half) / 6.086));
	    double thf = 1.0 / (1.432e-5 * exp(-(state[0] + 1.196) / 6.285) + 6.149 * exp((state[0] + 0.5096) / 20.27));
	    double ths = 1.0 / (0.009794 * exp(-(state[0] + 17.95) / 28.05) + 0.3343 * exp((state[0] + 5.730) / 56.66));
	    double Ahf = 0.99;
	    double Ahs = 1.0 - Ahf;
	    state[10] = hss - (hss - state[10]) * exp(-dt / thf);
	    state[11] = hss - (hss - state[11]) * exp(-dt / ths);
	    double h = Ahf * state[10] + Ahs * state[11];
	    double jss = hss;
	    double tj = 2.038 + 1.0 / (0.02136 * exp(-(state[0] + 100.6) / 8.281) + 0.3052 * exp((state[0] + 0.9941) / 38.45));
	    state[12] = jss - (jss - state[12]) * exp(-dt / tj);
	    double hssp = 1.0 / (1 + exp((state[0] + 89.1 + INa_CaMh_half) / 6.086));
	    double thsp = 3.0 * ths;
	    state[13] = hssp - (hssp - state[13]) * exp(-dt / thsp);
	    double hp = Ahf * state[10] + Ahs * state[13];
	    double tjp = 1.46 * tj;
	    state[14] = jss - (jss - state[14]) * exp(-dt / tjp);
	    // printf("debug _ 1\n");
	    double GNa = 75;
	    double fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
	    INa = GNa_ISO * GNa * (state[0] - ENa) * state[9] * state[9] * state[9] * ((1.0 - fINap) * h * state[12] + fINap * hp * state[14]);
	}*/



	/*using m gate of TNNP INa model in the ORd INa model*/
	double tm = 1.0 / (6.765 * exp((state[0] + 11.64) / 34.77) + 8.552 * exp(-(state[0] + 77.42) / 5.955));
	/*{
		double mss = 1. / ((1. + exp((-56.86 - state[0]) / 9.03)) * (1. + exp((-56.86 - state[0]) / 9.03)));
		state[9] = mss - (mss - state[9]) * exp(-dt / tm);
		double hss = 1.0 / (1 + exp((state[0] + 82.90 + INah_half) / 6.086));
		double thf = 1.0 / (1.432e-5 * exp(-(state[0] + 1.196) / 6.285) + 6.149 * exp((state[0] + 0.5096) / 20.27));
		double ths = 1.0 / (0.009794 * exp(-(state[0] + 17.95) / 28.05) + 0.3343 * exp((state[0] + 5.730) / 56.66));
		double Ahf = 0.99;
		double Ahs = 1.0 - Ahf;
		state[10] = hss - (hss - state[10]) * exp(-dt / thf);
		state[11] = hss - (hss - state[11]) * exp(-dt / ths);
		double h = Ahf * state[10] + Ahs * state[11];
		double jss = hss;
		double tj = 2.038 + 1.0 / (0.02136 * exp(-(state[0] + 100.6) / 8.281) + 0.3052 * exp((state[0] + 0.9941) / 38.45));
		state[12] = jss - (jss - state[12]) * exp(-dt / tj);
		double hssp = 1.0 / (1 + exp((state[0] + 89.1 + INa_CaMh_half) / 6.086));
		double thsp = 3.0 * ths;
		state[13] = hssp - (hssp - state[13]) * exp(-dt / thsp);
		double hp = Ahf * state[10] + Ahs * state[13];
		double tjp = 1.46 * tj;
		state[14] = jss - (jss - state[14]) * exp(-dt / tjp);
		// printf("debug _ 1\n");
		double GNa = 75;
		double fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
		INa = GNa_ISO * GNa * (state[0] - ENa) * state[9] * state[9] * state[9] * ((1.0 - fINap) * h * state[12] + fINap * hp * state[14]);
	}
	*/

	// Reversal Potentials
	// Eca = 13.35 * log(CRN_cac / cai);




	// //INa
	double GNa = 10.1;
	double V = state[0];

	double INa_BA = state[41];
	double INa_BI = state[42];
	double m = state[9];
	double h = state[10];
	double j = state[11];

	INa = 1 * GNa * (1 - INa_BA - INa_BI) * m * m * m * h * j * (V - ENa); // drug blocked INa
	het.INa = INa / 1;


	double dBAdt = het.drug_INa_Ka * het.drug_INa_concen * m * m * m * h * j * ( 1 - INa_BA - INa_BI) - het.drug_INa_La * INa_BA ; // * exp(-z*V*F/R/T)
	double dBIdt = het.drug_INa_Ki * het.drug_INa_concen * (1 - h) * (1 - INa_BA - INa_BI) - het.drug_INa_Li * INa_BI; // * exp(-z*V*F/R/T)


	INa_BA += dt * dBAdt;   // blockade to Open gate.
	INa_BI += dt * dBIdt;	// Blockade to closed states.
	state[41] = INa_BA;
	state[42] = INa_BI;

	// m gate
	double alpha =  0.32 * ((V) + 47.13) / (1 - exp(-0.1 * ((V) + 47.13)));
	if (fabs(V + 47.13) < 1e-10) alpha = 3.2;
	double beta = 0.08 * exp(-(V) / 11.0);

	double tau = 1 / (alpha + beta);
	double inf = alpha * tau;


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

	state[9] = m;
	state[10] = h;
	state[11] = j;


	/*TNNP INa formula*/

	/*double mss =1. / ((1. + exp((-56.86 - state[0]) / 9.03)) * (1. + exp((-56.86 - state[0]) / 9.03)));
	double tm = 1.0 / (6.765 * exp((state[0] + 11.64) / 34.77) + 8.552 * exp(-(state[0] + 77.42) / 5.955));
	state[9] = mss - (mss - state[9]) * exp(-dt / tm);


	double hss = 1. / ((1. + exp((state[0] + 71.55+ INah_half) / 7.43)) * (1. + exp((state[0] + 71.55+ INah_half) / 7.43)));
	double thf = 1.0 / (1.432e-5 * exp(-(state[0] + 1.196) / 6.285) + 6.149 * exp((state[0] + 0.5096) / 20.27));
	double ths = 1.0 / (0.009794 * exp(-(state[0] + 17.95) / 28.05) + 0.3343 * exp((state[0] + 5.730) / 56.66));
	double Ahf = 0.99;
	double Ahs = 1.0 - Ahf;
	state[10] = hss - (hss - state[10]) * exp(-dt / thf);
	state[11] = hss - (hss - state[11]) * exp(-dt / ths);
	double h = Ahf * state[10] + Ahs * state[11];
	double jss = hss;
	double tj = 2.038 + 1.0 / (0.02136 * exp(-(state[0] + 100.6) / 8.281) + 0.3052 * exp((state[0] + 0.9941) / 38.45));
	state[12] = jss - (jss - state[12]) * exp(-dt / tj);
	double hssp = 1.0 / (1 + exp((state[0] + 89.1 + INa_CaMh_half) / 6.086));
	double thsp = 3.0 * ths;
	state[13] = hssp - (hssp - state[13]) * exp(-dt / thsp);
	double hp = Ahf * state[10] + Ahs * state[13];
	double tjp = 1.46 * tj;
	state[14] = jss - (jss - state[14]) * exp(-dt / tjp);
	// printf("debug _ 1\n");
	double GNa = 75;
	double fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
	INa = GNa_ISO * GNa * (state[0] - ENa) * state[9] * state[9] * state[9] * ((1.0 - fINap) * h * state[12] + fINap * hp * state[14]);*/


	double tmL = tm;
	{	double mLss = 1.0 / (1.0 + exp((-(state[0] + 42.85)) / 5.264));
		state[15] = mLss - (mLss - state[15]) * exp(-dt / tmL);
		double hLss = 1.0 / (1.0 + exp((state[0] + 87.61) / 7.488));
		double thL = 200.0;
		state[16] = hLss - (hLss - state[16]) * exp(-dt / thL);
		double hLssp = 1.0 / (1.0 + exp((state[0] + 93.81) / 7.488));
		double thLp = 3.0 * thL;
		state[17] = hLssp - (hLssp - state[17]) * exp(-dt / thLp);
		double GNaL = 0.0075;
		if (celltype == LVEPI)
		{
			GNaL *= 0.6;
		}
		double fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
		INaL = GNaL * (state[0] - ENa) * state[15] * ((1.0 - fINaLp) * state[16] + fINaLp * state[17]);
	}

	{	double ass = 1.0 / (1.0 + exp((-(state[0] - 14.34)) / 14.82));
		double ta = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(state[0] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((state[0] + 100.0) / 29.3814)));
		state[18] = ass - (ass - state[18]) * exp(-dt / ta);
		double iss = 1.0 / (1.0 + exp((state[0] + 43.94 + InAcTo_Vhalf_ABh) / 5.711));
		double delta_LVEPI;
		if (celltype == LVEPI)
		{
			delta_LVEPI = 1.0 - (0.95 / (1.0 + exp((state[0] + 70.0) / 5.0)));
		}
		else
		{
			delta_LVEPI = 1.0;
		}
		double tiF = 4.562 + 1 / (0.3933 * exp((-(state[0] + 100.0)) / 100.0) + 0.08004 * exp((state[0] + 50.0) / 16.59));
		double tiS = 23.62 + 1 / (0.001416 * exp((-(state[0] + 96.52)) / 59.05) + 1.780e-8 * exp((state[0] + 114.1) / 8.079));
		tiF *= delta_LVEPI;
		tiS *= delta_LVEPI;
		double AiF = 1.0 / (1.0 + exp((state[0] - 213.6) / 151.2));
		double AiS = 1.0 - AiF;
		state[19] = iss - (iss - state[19]) * exp(-dt / tiF);
		state[20] = iss - (iss - state[20]) * exp(-dt / tiS);
		double i = AiF * state[19] + AiS * state[20];
		double assp = 1.0 / (1.0 + exp((-(state[0] - 24.34)) / 14.82));
		state[21] = assp - (assp - state[21]) * exp(-dt / ta);
		double dti_develop = 1.354 + 1.0e-4 / (exp((state[0] - 167.4) / 15.89) + exp(-(state[0] - 12.23) / 0.2154));
		double dti_recover = 1.0 - 0.5 / (1.0 + exp((state[0] + 70.0) / 20.0));
		double tiFp = dti_develop * dti_recover * tiF;
		double tiSp = dti_develop * dti_recover * tiS;
		state[22] = iss - (iss - state[22]) * exp(-dt / tiFp);
		state[23] = iss - (iss - state[23]) * exp(-dt / tiSp);
		double ip = AiF * state[22] + AiS * state[23];
		double Gto = 0.02;
		if (celltype == LVEPI)
		{
			Gto *= 4.0;
			// printf("LVEPI\n");
		}
		if (celltype == LVMCELL)
		{
			Gto *= 4.0;
			// printf("LVMCELL\n");
		}
		double fItop = (1.0 / (1.0 + KmCaMK / CaMKa));
		Ito = GTo_ABh * Gto * (state[0] - EK) * ((1.0 - fItop) * state[18] * i + fItop * state[21] * ip);
	}
	// printf("debug _ 2\n");

	{	double dss = 1.0 / (1.0 + exp((-(state[0] + 3.940 + ICaLd_half)) / 4.230));


		double td = 0.6 + 1.0 / (exp(-0.05 * (state[0] + 6.0)) + exp(0.09 * (state[0] + 14.0)));
		state[24] = dss - (dss - state[24]) * exp(-dt / td);

		// double fss = 1.0 / (1.0 + exp((state[0] + 19.58 + ICaLf_half) / 3.696));  // ref 1
		// double fss = 1.0 / (1.0 + exp((state[0] + 20.5 + ICaLf_half) / 3.1775));  // ref 1
		double fss = 1.0 / (1.0 + exp((state[0] + 20.25 + ICaLf_half) / 3.696));  // ref 1

		// double fss = 1.0 / (1.0 + exp((state[0] + 19.58 + ICaLf_half) / 3.1775));  // ref 1
		// double fss = (1.0 / (1.0 + exp((state[0] + 21.979 + ICaLf_half) / 3.1775)) + 0.0001) * (1.0 / 1.0001); // ref 2
		// double fss = (1.0 / (1.0 + exp((state[0] + 20.979 + ICaLf_half) / 3.1775)) ) ; // ref 2




		double tff = Tau_ICaL_ff_ISO * (7.0 + 1.0 / (0.0045 * exp(-(state[0] + 20.0) / 10.0) + 0.0045 * exp((state[0] + 20.0) / 10.0)));
		double tfs = 1000.0 + 1.0 / (0.000035 * exp(-(state[0] + 5.0) / 4.0) + 0.000035 * exp((state[0] + 5.0) / 6.0));
		double Aff = 0.6;
		double Afs = 1.0 - Aff;
		state[25] = fss - (fss - state[25]) * exp(-dt / tff);
		state[26] = fss - (fss - state[26]) * exp(-dt / tfs);
		double f = Aff * state[25] + Afs * state[26];
		double fcass = fss;
		double tfcaf = 7.0 + 1.0 / (0.04 * exp(-(state[0] - 4.0) / 7.0) + 0.04 * exp((state[0] - 4.0) / 7.0));
		double tfcas = 100.0 + 1.0 / (0.00012 * exp(-state[0] / 3.0) + 0.00012 * exp(state[0] / 7.0));
		double Afcaf = 0.3 + 0.6 / (1.0 + exp((state[0] - 10.0) / 10.0));
		double Afcas = 1.0 - Afcaf;
		state[27] = fcass - (fcass - state[27]) * exp(-dt / tfcaf);
		state[28] = fcass - (fcass - state[28]) * exp(-dt / tfcas);
		double fca = Afcaf * state[27] + Afcas * state[28];
		double tjca = 75.0;
		state[29] = fcass - (fcass - state[29]) * exp(-dt / tjca);
		double tffp = 2.5 * tff;
		state[31] = fss - (fss - state[31]) * exp(-dt / tffp);
		double fp = Aff * state[31] + Afs * state[26];
		double tfcafp = 2.5 * tfcaf;
		state[32] = fcass - (fcass - state[32]) * exp(-dt / tfcafp);
		double fcap = Afcaf * state[32] + Afcas * state[28];
		double Kmn = 0.002;
		double k2n = 1000.0;
		double km2n = state[29] * 1.0;
		double anca = 1.0 / (k2n / km2n + pow(1.0 + Kmn / state[6], 4.0));
		state[30] = anca * k2n / km2n - (anca * k2n / km2n - state[30]) * exp(-km2n * dt);

		double Cass = state[6];

		if (Cass > 0.03)
		{
			Cass = 0.03;   // 0.03mM ceiling, from ref 1. * 1. Ohara et al. Arrhythmia formation in subclinical (“silent”) long QT syndrome requires multiple insults: Quantitative mechanistic study using the KCNQ1 mutation Q357R as example
		}
		double PhiCaL = 4.0 * vffrt * (Cass * exp(2.0 * vfrt) - 0.341 * cao) / (exp(2.0 * vfrt) - 1.0);
		double PhiCaNa = 1.0 * vffrt * (0.75 * state[2] * exp(1.0 * vfrt) - 0.75 * nao) / (exp(1.0 * vfrt) - 1.0);
		double PhiCaK = 1.0 * vffrt * (0.75 * state[4] * exp(1.0 * vfrt) - 0.75 * ko) / (exp(1.0 * vfrt) - 1.0);

		double PCa = Pca_ISO * 0.0001;
		if (celltype == LVEPI)
		{
			PCa *= 1.2;
		}
		if (celltype == LVMCELL)
		{
			PCa *= 2.5;
		}
		double PCap = 1.1 * PCa;
		double PCaNa = 0.00125 * PCa;
		double PCaK = 3.574e-4 * PCa;
		double PCaNap = 0.00125 * PCap;
		double PCaKp = 3.574e-4 * PCap;
		double fICaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
		ICaL = (1.0 - fICaLp) * PCa * PhiCaL * state[24] * (f * (1.0 - state[30]) + state[29] * fca * state[30]) + fICaLp * PCap * PhiCaL * state[24] * (fp * (1.0 - state[30]) + state[29] * fcap * state[30]);
		ICaNa = (1.0 - fICaLp) * PCaNa * PhiCaNa * state[24] * (f * (1.0 - state[30]) + state[29] * fca * state[30]) + fICaLp * PCaNap * PhiCaNa * state[24] * (fp * (1.0 - state[30]) + state[29] * fcap * state[30]);
		ICaK = (1.0 - fICaLp) * PCaK * PhiCaK * state[24] * (f * (1.0 - state[30]) + state[29] * fca * state[30]) + fICaLp * PCaKp * PhiCaK * state[24] * (fp * (1.0 - state[30]) + state[29] * fcap * state[30]);
	}
	{	double xrss = 1.0 / (1.0 + exp((-(state[0] + 8.337)) / 6.789));
		double txrf = 12.98 + 1.0 / (0.3652 * exp((state[0] - 31.66) / 3.869) + 4.123e-5 * exp((-(state[0] - 47.78)) / 20.38));
		double txrs = 1.865 + 1.0 / (0.06629 * exp((state[0] - 34.70) / 7.355) + 1.128e-5 * exp((-(state[0] - 29.74)) / 25.94));
		double Axrf = 1.0 / (1.0 + exp((state[0] + 54.81) / 38.21));
		double Axrs = 1.0 - Axrf;
		state[33] = xrss - (xrss - state[33]) * exp(-dt / txrf);
		state[34] = xrss - (xrss - state[34]) * exp(-dt / txrs);
		double xr = Axrf * state[33] + Axrs * state[34];
		double rkr = 1.0 / (1.0 + exp((state[0] + 55.0) / 75.0)) * 1.0 / (1.0 + exp((state[0] - 10.0) / 30.0));
		double GKr = 0.046;
		if (celltype == LVEPI)
		{
			GKr *= 1.3;
		}
		if (celltype == LVMCELL)
		{
			// GKr *= 0.8;
			GKr *= 0.8;  // 1.3*0.91 in later paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
		}
		IKr = GKr * sqrt(ko / 5.4) * xr * rkr * (state[0] - EK);
	}


	if (mutation == ORd)
	{
		double xs1ss = 1.0 / (1.0 + exp((-(state[0] + IKs_ac_shift_ISO + 11.60)) / 8.932));
		double txs1 = Tau_Ks * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;
		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.0193 * exp((-(state[0] + 66.54)) / 31.0));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034;
		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);
	} else if (mutation == IKS_WT)
	{

		double xs1ss = 1.0 / (1.0 + exp((-(state[0] + 2.4 + IKs_ac_shift_ISO)) / 13.2));
		// in the new IKs model, txs1 was timed by 0.5 to better reproduce the results of voltage clamp.

		double txs1 = Tau_Ks * 0.5 * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;


		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.0071038 * exp((-(state[0] + 90.5741763432)) / 69.080));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034 * 1.2; // to make sure that the new IKS_WT IKs produced same APD90 as the original ORd IKs model
		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);
	} else if (mutation == HET)
	{
		double xs1ss = 1.0 / (1.0 + exp((-(state[0] - 25.4 + IKs_ac_shift_ISO)) / 17.8));
		double txs1 = Tau_Ks * 0.5 * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;


		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.017671 * exp((-(state[0] + 93.997 )) / 43.8948649349));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034 * 1.2; // to make sure that the new IKS_WT IKs produced same APD90 as the original ORd IKs model,
		// and the conductance of IKs has been scaled so that the tail current at 30mV in the voltage clamp matches ratio HET/IKS_WT in the original experimetnal data
		// see Figure 2E in paper Wu et al. 2014. A Molecular Mechanism for Adrenergic-Induced Long QT Syndrome  Jie Wu, P H D,*yz Nobu Naiki, MD,y Wei-Guang Ding, MD, P H D,z

		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);

	} else if (mutation == G269S)
	{
		double xs1ss = 1.0 / (1.0 + exp((-(state[0] - 68.3 + IKs_ac_shift_ISO)) / 23.4));
		double txs1 = Tau_Ks * 0.5 * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;

		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.0290693 * exp((-(state[0] +  106.263263059)) / 48.7397674066));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034 * 5.2 * 1.2;
		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);
	} else {
		printf("wrong IKs mutation option!!!\n");
		exit(0);
	}



	/* effects of Acacetin */
	
	IKs =IKs*het.GKs;
	Ito=Ito*het.Gto;
	IKr = IKr*het.GKr;



	double xk1ss = 1.0 / (1.0 + exp(-(state[0] + 2.5538 * ko + 144.59) / (1.5692 * ko + 3.8115)));
	double txk1 = 122.2 / (exp((-(state[0] + 127.2)) / 20.36) + exp((state[0] + 236.8) / 69.33));
	state[37] = xk1ss - (xk1ss - state[37]) * exp(-dt / txk1);
	double rk1 = 1.0 / (1.0 + exp((state[0] + 105.8 - 2.6 * ko) / 9.493));
	double GK1 = 0.1908;
	if (celltype == LVEPI)
	{
		GK1 *= 1.2;
	}
	if (celltype == LVMCELL)
	{
		GK1 *= 1.3;
	}
	IK1 = GK1 * sqrt(ko) * rk1 * state[37] * (state[0] - EK);

	// double 15.0 = 15.0;
	// double 5.0 = 5.0;
	// double 88.12 = 88.12;
	// double 12.5 = 12.5;
	double x1;
	double x2;
	double x3;
	double x4;
	double E1, E2, E3, E4;
	{	double wna = 6.0e4;
		double wca = 6.0e4;
		double wnaca = 5.0e3;
		double kcaon = 1.5e6;
		double kcaoff = 5.0e3;
		double qna = 0.5224;
		double qca = 0.1670;
		double hca = exp((qca * state[0] * F) / (R * T));
		double hna = exp((qna * state[0] * F) / (R * T));
		double h1 = 1 + state[1] / 88.12 * (1 + hna);
		double h2 = (state[1] * hna) / (88.12 * h1);
		double h3 = 1.0 / h1;
		double h4 = 1.0 + state[1] / 15.0 * (1 + state[1] / 5.0);
		double h5 = state[1] * state[1] / (h4 * 15.0 * 5.0);
		double h6 = 1.0 / h4;
		double h7 = 1.0 + nao / 88.12 * (1.0 + 1.0 / hna);
		double h8 = nao / (88.12 * hna * h7);
		double h9 = 1.0 / h7;
		double h10 = 12.5 + 1.0 + nao / 15.0 * (1.0 + nao / 5.0);
		double h11 = nao * nao / (h10 * 15.0 * 5.0);
		double h12 = 1.0 / h10;
		double k1 = h12 * cao * kcaon;
		double k2 = kcaoff;
		double k3p = h9 * wca;
		double k3pp = h8 * wnaca;
		double k3 = k3p + k3pp;
		double k4p = h3 * wca / hca;
		double k4pp = h2 * wnaca;
		double k4 = k4p + k4pp;
		double k5 = kcaoff;
		double k6 = h6 * state[5] * kcaon;
		double k7 = h5 * h2 * wna;
		double k8 = h8 * h11 * wna;
		x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
		x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
		x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
		x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
		E1 = x1 / (x1 + x2 + x3 + x4);
		E2 = x2 / (x1 + x2 + x3 + x4);
		E3 = x3 / (x1 + x2 + x3 + x4);
		E4 = x4 / (x1 + x2 + x3 + x4);
		double KmCaAct = 150.0e-6;
		double allo = 1.0 / (1.0 + pow(KmCaAct / state[5], 2.0));

		double JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
		double JncxCa = E2 * k2 - E1 * k1;
		double Gncx = 0.0008;
		if (celltype == LVEPI)
		{
			Gncx *= 1.1;
		}
		if (celltype == LVMCELL)
		{
			Gncx *= 1.4;
		}
		INaCa_i = 0.8 * Gncx * allo * (zna * JncxNa + zca * JncxCa);

		h1 = 1 + state[2] / 88.12 * (1 + hna);
		h2 = (state[2] * hna) / (88.12 * h1);
		h3 = 1.0 / h1;
		h4 = 1.0 + state[2] / 15.0 * (1 + state[2] / 5.0);
		h5 = state[2] * state[2] / (h4 * 15.0 * 5.0);
		h6 = 1.0 / h4;
		h7 = 1.0 + nao / 88.12 * (1.0 + 1.0 / hna);
		h8 = nao / (88.12 * hna * h7);
		h9 = 1.0 / h7;
		h10 = 12.5 + 1.0 + nao / 15.0 * (1 + nao / 5.0);
		h11 = nao * nao / (h10 * 15.0 * 5.0);
		h12 = 1.0 / h10;
		k1 = h12 * cao * kcaon;
		k2 = kcaoff;
		k3p = h9 * wca;
		k3pp = h8 * wnaca;
		k3 = k3p + k3pp;
		k4p = h3 * wca / hca;
		k4pp = h2 * wnaca;
		k4 = k4p + k4pp;
		k5 = kcaoff;
		k6 = h6 * state[6] * kcaon;
		k7 = h5 * h2 * wna;
		k8 = h8 * h11 * wna;
		x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
		x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
		x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
		x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
		E1 = x1 / (x1 + x2 + x3 + x4);
		E2 = x2 / (x1 + x2 + x3 + x4);
		E3 = x3 / (x1 + x2 + x3 + x4);
		E4 = x4 / (x1 + x2 + x3 + x4);
		KmCaAct = 150.0e-6;
		allo = 1.0 / (1.0 + pow(KmCaAct / state[6], 2.0));
		JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
		JncxCa = E2 * k2 - E1 * k1;
		INaCa_ss = 0.2 * Gncx * allo * (zna * JncxNa + zca * JncxCa);

		INaCa = INaCa_i + INaCa_ss;
	}
	// printf("debug _ 3\n");

	{	double k1p = 949.5;
		double k1m = 182.4;
		double k2p = 687.2;
		double k2m = 39.4;
		double k3p = 1899.0;
		double k3m = 79300.0;
		double k4p = 639.0;
		double k4m = 40.0;
		// double Knai0 = 9.073;
		// double Knao0 = 27.78;
		double delta = -0.1550;
		double Knai = NaK_ISO * (9.073 * exp((delta * state[0] * F) / (3.0 * R * T)));
		double Knao = 27.78 * exp(((1.0 - delta) * state[0] * F) / (3.0 * R * T));
		double Kki = 0.5;
		double Kko = 0.3582;
		double MgADP = 0.05;
		double MgATP = 9.8;
		double Kmgatp = 1.698e-7;
		double H = 1.0e-7;
		double eP = 4.2;
		double Khp = 1.698e-7;
		double Knap = 224.0;
		double Kxkur = 292.0;
		double P = eP / (1.0 + H / Khp + state[1] / Knap + state[3] / Kxkur);
		double a1 = (k1p * pow(state[1] / Knai, 3.0)) / (pow(1.0 + state[1] / Knai, 3.0) + pow(1.0 + state[3] / Kki, 2.0) - 1.0);
		double b1 = k1m * MgADP;
		double a2 = k2p;
		double b2 = (k2m * pow(nao / Knao, 3.0)) / (pow(1.0 + nao / Knao, 3.0) + pow(1.0 + ko / Kko, 2.0) - 1.0);
		double a3 = (k3p * pow(ko / Kko, 2.0)) / (pow(1.0 + nao / Knao, 3.0) + pow(1.0 + ko / Kko, 2.0) - 1.0);
		double b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp);
		double a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp);
		double b4 = (k4m * pow(state[3] / Kki, 2.0)) / (pow(1.0 + state[1] / Knai, 3.0) + pow(1.0 + state[3] / Kki, 2.0) - 1.0);
		double x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
		double x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
		double x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
		double x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
		E1 = x1 / (x1 + x2 + x3 + x4);
		E2 = x2 / (x1 + x2 + x3 + x4);
		E3 = x3 / (x1 + x2 + x3 + x4);
		E4 = x4 / (x1 + x2 + x3 + x4);
		double zk = 1.0;
		double JnakNa = 3.0 * (E1 * a3 - E2 * b3);
		double JnakK = 2.0 * (E4 * b1 - E3 * a1);
		double Pnak = 30;
		if (celltype == LVEPI)
		{
			Pnak *= 0.9;
		}
		if (celltype == LVMCELL)
		{
			Pnak *= 0.7;
		}
		INaK = Pnak * (zna * JnakNa + zk * JnakK);

		double xkb = 1.0 / (1.0 + exp(-(state[0] - 14.48) / 18.34));
		double GKb = 0.003;
		IKb = GKb_ISO * GKb * xkb * (state[0] - EK);

		double PNab = 3.75e-10;
		INab = PNab * vffrt * (state[1] * exp(vfrt) - nao) / (exp(vfrt) - 1.0);

		double PCab = 2.5e-8;
		ICab = PCab * 4.0 * vffrt * (state[5] * exp(2.0 * vfrt) - 0.341 * cao) / (exp(2.0 * vfrt) - 1.0);

		/*double GpCa = 0.005;
		IpCa = GpCa * state[5] / (0.00025 + state[5]);  // note the cai pump was increased to expel more cai.*/
		double GpCa = 0.0005;
		IpCa = GpCa * state[5] / (0.0005 + state[5]);
	}


	// // original..
	//  double GpCa = 0.0005;
	// IpCa = GpCa * cai / (0.0005 + cai);

	JdiffNa = (state[2] - state[1]) / 2.0;
	JdiffK = (state[4] - state[3]) / 2.0;
	Jdiff = (state[6] - state[5]) / 0.2;

	double bt = 4.75;
	double a_rel = A_rel_ISO  * 0.5 * bt ;

	double Jrel_inf = a_rel * (-ICaL) / (1.0 + pow(1.5 / state[8], 8.0));
	if (celltype == LVMCELL)
	{
		Jrel_inf *= 1.7;
	}
	double tau_rel = tau_rel_ISO * bt / (1.0 + 0.0123 / state[8]);
	if (tau_rel < 0.005)
	{
		tau_rel = 0.005;
	}
	state[38] = Jrel_inf - (Jrel_inf - state[38]) * exp(-dt / tau_rel);
	double btp = 1.25 * bt;
	double a_relp = A_rel_ISO * 0.5 * btp;
	double Jrel_infp = a_relp * (-ICaL) / (1.0 + pow(1.5 / state[8], 8.0));
	if (celltype == LVMCELL)
	{
		Jrel_infp *= 1.7;
	}
	double tau_relp = tau_rel_ISO * btp / (1.0 + 0.0123 / state[8]);
	if (tau_relp < 0.005)
	{
		tau_relp = 0.005;
	}
	state[39] = Jrel_infp - (Jrel_infp - state[39]) * exp(-dt / tau_relp);
	double fJrelp = (1.0 / (1.0 + KmCaMK / CaMKa));
	Jrel = (1.0 - fJrelp) * state[38] + fJrelp * state[39];




	double Jupnp = 0.004375 * state[5] / (state[5] + 0.00092 * Jup_ISO);
	double Jupp = 2.75 * 0.004375 * state[5] / (state[5] + ( 0.00092 - 0.00017) * Jup_ISO);


	if (celltype == LVEPI)
	{
		Jupnp *= 1.3;
		Jupp *= 1.3;
	}


	double fJupp = (1.0 / (1.0 + KmCaMK / CaMKa));
	Jleak = 0.0039375 * state[7] / 15.0;
	Jup = (1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak;

	Jtr = (state[7] - state[8]) / 100.0;

	state[1] += dt * (-(INa + INaL + 3.0 * INaCa_i + 3.0 * INaK + INab) * Acap / (F * vmyo) + JdiffNa * vss / vmyo);
	state[2] += dt * (-(ICaNa + 3.0 * INaCa_ss) * Acap / (F * vss) - JdiffNa);

	state[3] += dt * (-(Ito + IKr + IKs + IK1 + IKb + Ist - 2.0 * INaK) * Acap / (F * vmyo) + JdiffK * vss / vmyo);
	state[4] += dt * (-(ICaK) * Acap / (F * vss) - JdiffK);
	double kmtrpn_cell = kmtrpn * kmtrpn_ISO;
	double Bcai;
	if (celltype == LVEPI)
	{
		Bcai = 1.0 / (1.0 + 1.3 * cmdnmax * kmcmdn / pow(kmcmdn + state[5], 2.0) + trpnmax * kmtrpn_cell / pow(kmtrpn_cell + state[5], 2.0));
	}
	else
	{
		Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / pow(kmcmdn + state[5], 2.0) + trpnmax * kmtrpn_cell / pow(kmtrpn_cell + state[5], 2.0));
	}
	state[5] += dt * (Bcai * (-(IpCa + ICab - 2.0 * INaCa_i) * Acap / (2.0 * F * vmyo) - Jup * vnsr / vmyo + Jdiff * vss / vmyo));

	double Bcass = 1.0 / (1.0 + BSRmax * KmBSR / pow(KmBSR + state[6], 2.0) + BSLmax * KmBSL / pow(KmBSL + state[6], 2.0));
	state[6] += dt * (Bcass * (-(ICaL - 2.0 * INaCa_ss) * Acap / (2.0 * F * vss) + Jrel * vjsr / vss - Jdiff));

	state[7] += dt * (Jup - Jtr * vjsr / vnsr);

	double Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / pow(kmcsqn + state[8], 2.0));
	state[8] += dt * (Bcajsr * (Jtr - Jrel));

	return  (INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + IKb + IpCa + ICab + Ist);
	// f// printf(stderr, "Ina2=%f\n", INa);

	// double I = INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + IKb + IpCa + ICab + Ist;
	// f// printf(stderr, "I=%f\n", I);
	// printf("debug _ 4\n");

	// ion -> IKr = IKr;
	// ion ->  IKs = IKs;
	// ion -> IK1 = IK1;
	// ion -> Ito = Ito;
	// // ion -> Isus = Isus;
	// ion -> INa = INa;
	// ion ->  IbNa = INab;
	// ion ->ICaL =  ICaL;
	// ion -> IbCa = ICab;
	// ion -> INaK = INaK;
	// ion ->  INaCa = INaCa;//!!!!!
	// ion -> IpCa = IpCa;
	// ion -> Cai = state[5];//!!!!

	// ion -> IpK = IpK;
	// ion -> If = If;
	// printf("debug _ 5\n");


	/*static double t;
	t += dt;
	static int count;
	count += 1;
	FILE * output = fopen("current.txt", "a");
	if (count % 10 == 0)
	{

	    fprintf(output, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t,Jup, Jrel, INa , INaL , Ito , ICaL , ICaNa , ICaK , IKr , IKs , IK1 , INaCa , INaK , INab , IKb , IpCa , ICab , Ist);
	}
	fclose(output);*/

	// IKs_rec = IKs;


}


double ORd_Model(double dt, double stims, double *state, TypeCell celltype, int ISO, IKsMutation mutation, float ABIndex)
{
	double ENa, EK, EKs;
	double INa, INaL, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCa_i, INaCa_ss, INaCa, INaK, IKb, INab, IpCa, ICab, Ist;
	double Jrel, Jup, Jtr, Jdiff, JdiffNa, JdiffK, Jleak;
	double CaMKa, CaMKb;
	Ist = stims;

	double INah_half, INa_CaMh_half, GNa_ISO;

	double ICaLd_half, ICaLf_half, Pca_ISO;
	double Tau_ICaL_ff_ISO;

	double Tau_Ks, GKs_ISO, IKs_ac_shift_ISO;

	double NaK_ISO;

	double GKb_ISO;


	double A_rel_ISO, tau_rel_ISO;
	double Jup_ISO;
	double kmtrpn_ISO;

	INah_half = INa_CaMh_half = 0.0;
	GNa_ISO = 1.0;
	ICaLd_half = 0.0;
	ICaLf_half = 0.0;
	Tau_ICaL_ff_ISO = 1.0;
	Pca_ISO = 1.0;
	GKs_ISO = 1.0;
	Tau_Ks = 1.0;
	NaK_ISO = 1.0;//
	GKb_ISO = 1.0;
	A_rel_ISO = 1.0;   // ref 2
	tau_rel_ISO = 1.0;  // ref 2
	Jup_ISO = 1.0;
	kmtrpn_ISO = 1.0;

	CaMKb = CaMKo * (1.0 - state[40]) / (1.0 + KmCaM / state[6]);
	CaMKa = CaMKb + state[40];
	state[40] += dt * (aCaMK * CaMKb * (CaMKb + state[40]) - bCaMK * state[40]);




	if (ISO == 1)
	{
		double time_1;
		/*        if (t > 100000)
		        {
		            time_1 = t - 100000;

		        }
		        else {
		            time_1 = 0.0;
		        }
		*/
		double value_1 = 1.0;// - exp(-time_1 / 7700);
		double value_2 = 1.0;// - exp(-time_1 / 39700);

		/*if (time_1 > 100000)
		{
		    value_2 = value_1 = 1.0;
		}*/
		INah_half = INa_CaMh_half = value_1 * 5.0;
		GNa_ISO =   1.0 + value_1 * 1.70;
		ICaLd_half = value_1 * 9.9; // ref 1 and 2 16.0
		// change to 8 will give good results  // the original value from rudy paper was 16.0
		/*Figure 3 in paper         // from  paper : Chen et al. L-Type Ca2+ Channel Density and Regulation Are Altered in Failing Human Ventricular Myocytes and Recover After Support With Mechanical Assist Devices
		ISO caused a significant leftward shift of the voltage dependence of I Ca,L activation in
		NF- (change of V 0.5 (dρ): -9.9 mV; P<0.05) and LVAD- (change of V 0.5 (dρ): -5.8 mV; P0.05) HVMs but no significant shift in F- (change of V 0.5 (dρ): -2.2 mV; PϾ0.05) HVMs.*/
		// td    = 0.6 + 1.0 / (exp(-0.05 * (v + 6.0 + 16.0)) + exp(0.09 * (v + 14.0 + 16.0)));  // in the rudy paper, the td was not changed by pka..
		Tau_ICaL_ff_ISO = 1 - value_1 * 0.5; // ISO accerartes the fast inactivation from paper chen et al. L-Type Ca2+ Channel Density and Regulation Are Altered in Failing Human Ventricular Myocytes and Recover After Support With Mechanical Assis
		/* In the paper */
		ICaLf_half = value_1 * 8.0; //ref 1
		Pca_ISO =  1.0 + value_1 * 1.5;

		GKs_ISO = 1.0 + value_2 * 2.2;
		IKs_ac_shift_ISO = value_2 * 7.0; // Impaired IKs channel activation by Ca2 +-dependent PKC shows correlation with emotion/arousal-triggered events in LQT1

		if (mutation == ORd or mutation == IKS_WT)
		{
			GKs_ISO = 1.0 + value_2 * 2.2;
			Tau_Ks = 1.0 - value_2 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
			IKs_ac_shift_ISO = 7.0 * value_2;
		} else if (mutation == HET)
		{
			GKs_ISO = 1 + value_2 * 0.88;
			Tau_Ks = 1.0 - value_2 * 0.4 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
			IKs_ac_shift_ISO = 7.0 * 0.4 * value_2;
			// GKs_ISO = 3.2 * 0.4;
		} else {
			// GKs_ISO = 3.2 * 0.4;
			// GKs_ISO = 2.23;
			// GKs_ISO = 1.0;
			GKs_ISO = 1 + value_2 * 0.88;
			Tau_Ks = 1.0 - value_2 * 0.4 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
			IKs_ac_shift_ISO = 7.0 * 0.4 * value_2;

		}
		// / say: previous work in human [(20) using original current formulations (21) with τxs1,PKA = 0.6 * τxs1 and GKs,PKA = 3.2 * GKs,
		NaK_ISO = 1 - value_1 * 0.3; //
		GKb_ISO = 1 + value_1 * 2.62; // ref 2 3.62 and ref 1 2.5
		A_rel_ISO = 1 + value_1 * 0.8; // ref 2
		tau_rel_ISO = 1 - value_1 * 0.6; // ref 2
		Jup_ISO = 1 - value_1 * 0.46;
		kmtrpn_ISO = 1.0 + value_1 * 0.6;
		/*GKs_ISO = 1.0 + value_2 * 2.2;
		Tau_Ks = 1.0 - value_2 * 0.4; // 1/1.2 -> 0.6    // in the paper says: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
		IKs_ac_shift_ISO = 7.0 * value_2;*/
	}


	/*  default values for ABIndex = 0.0 (Value set 0 == Apical Cells)
	        and ABIndex ranges from 0.0 -> 1.0 */
	/*
	*   0.0 -- Apical Cells
	*   1.0 -- Basal Cells
	*/

	double GTo_ABh = 1.0;
	double InAcTo_Vhalf_ABh = 0.0;
	double InAcTo_Vk_ABh = 1.0;
	double GKs_ABh = 1.0;
	double TauKs_ABh = 1.0;
	double AcKs_Vhalf_ABh = 0.0;
	double AB_IKr_ActVhalf = 0.0;




	if (ABIndex < 0.0 or ABIndex > 1.0) {
		std::cerr << "Wrong Apcial-Basal Ratio, at ABIndex = " << ABIndex << std:: endl;
		std::cerr << "Program Existing... " << std::endl;
		std::exit(0);
	}
	/* V1, ABIndex   == 0.0 -> Apical Cells */
	GTo_ABh          =  1.0 - ABIndex * (1 - (16.5) / 29.6) ; //(1.0 + (1 - ABIndex) * (/ 16.529.6  - 1));      // Apical cells: 1.0; Basal Cells: 16.5/29.6
	InAcTo_Vhalf_ABh = ABIndex * 4.0;                          // Apical Cells: 0; Basal Cells: 4.0;
	// InAcTo_Vk_ABh    = 1.0 + (1 - ABIndex) * ((4.5 - 3.4) / 3.4);   // Apical cells: 1.0; Basal Cells: 1.0  // not used here..
	GKs_ABh          = 1 - ABIndex * (1.0 - 2.1 / 5.6);  // Apical Cells: 1.0; Basal Cells: 2.1/5.6
	// TauKs_ABh        = 1.0 + (1 - ABIndex) * ((358.0 - 516.0) / 516.0);      //  // not used here.. refer to TNNP model

	// Iion_out *ion;
	// printf("%d\n",celltype);

	// double dt =  dti;
	// revpots(state);
	ENa = (R * T / F) * log(nao / state[1]);
	EK = (R * T / F) * log(ko / state[3]);
	EKs = (R * T / F) * log((ko + 0.01833 * nao) / (state[3] + 0.01833 * state[1]));

	// f// printf(stderr, "Ena=%f state10=%f\n", ENa, state[10]);
	// RGC(dt, state, celltype);
	CaMKb = CaMKo * (1.0 - state[40]) / (1.0 + KmCaM / state[6]);
	CaMKa = CaMKb + state[40];
	double vffrt = state[0] * F * F / (R * T);
	double vfrt = state[0] * F / (R * T);


	/* ORd INa*/

	/*double tm = 1.0 / (6.765 * exp((state[0] + 11.64) / 34.77) + 8.552 * exp(-(state[0] + 77.42) / 5.955));
	{   double mss = 1.0 / (1.0 + exp((-(state[0] + 39.57)) / 9.871));
	    state[9] = mss - (mss - state[9]) * exp(-dt / tm);
	    double hss = 1.0 / (1 + exp((state[0] + 82.90 + INah_half) / 6.086));
	    double thf = 1.0 / (1.432e-5 * exp(-(state[0] + 1.196) / 6.285) + 6.149 * exp((state[0] + 0.5096) / 20.27));
	    double ths = 1.0 / (0.009794 * exp(-(state[0] + 17.95) / 28.05) + 0.3343 * exp((state[0] + 5.730) / 56.66));
	    double Ahf = 0.99;
	    double Ahs = 1.0 - Ahf;
	    state[10] = hss - (hss - state[10]) * exp(-dt / thf);
	    state[11] = hss - (hss - state[11]) * exp(-dt / ths);
	    double h = Ahf * state[10] + Ahs * state[11];
	    double jss = hss;
	    double tj = 2.038 + 1.0 / (0.02136 * exp(-(state[0] + 100.6) / 8.281) + 0.3052 * exp((state[0] + 0.9941) / 38.45));
	    state[12] = jss - (jss - state[12]) * exp(-dt / tj);
	    double hssp = 1.0 / (1 + exp((state[0] + 89.1 + INa_CaMh_half) / 6.086));
	    double thsp = 3.0 * ths;
	    state[13] = hssp - (hssp - state[13]) * exp(-dt / thsp);
	    double hp = Ahf * state[10] + Ahs * state[13];
	    double tjp = 1.46 * tj;
	    state[14] = jss - (jss - state[14]) * exp(-dt / tjp);
	    // printf("debug _ 1\n");
	    double GNa = 75;
	    double fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
	    INa = GNa_ISO * GNa * (state[0] - ENa) * state[9] * state[9] * state[9] * ((1.0 - fINap) * h * state[12] + fINap * hp * state[14]);
	}*/



	/*using m gate of TNNP INa model in the ORd INa model*/
	double tm = 1.0 / (6.765 * exp((state[0] + 11.64) / 34.77) + 8.552 * exp(-(state[0] + 77.42) / 5.955));
	{
		double mss = 1. / ((1. + exp((-56.86 - state[0]) / 9.03)) * (1. + exp((-56.86 - state[0]) / 9.03)));
		state[9] = mss - (mss - state[9]) * exp(-dt / tm);
		double hss = 1.0 / (1 + exp((state[0] + 82.90 + INah_half) / 6.086));
		double thf = 1.0 / (1.432e-5 * exp(-(state[0] + 1.196) / 6.285) + 6.149 * exp((state[0] + 0.5096) / 20.27));
		double ths = 1.0 / (0.009794 * exp(-(state[0] + 17.95) / 28.05) + 0.3343 * exp((state[0] + 5.730) / 56.66));
		double Ahf = 0.99;
		double Ahs = 1.0 - Ahf;
		state[10] = hss - (hss - state[10]) * exp(-dt / thf);
		state[11] = hss - (hss - state[11]) * exp(-dt / ths);
		double h = Ahf * state[10] + Ahs * state[11];
		double jss = hss;
		double tj = 2.038 + 1.0 / (0.02136 * exp(-(state[0] + 100.6) / 8.281) + 0.3052 * exp((state[0] + 0.9941) / 38.45));
		state[12] = jss - (jss - state[12]) * exp(-dt / tj);
		double hssp = 1.0 / (1 + exp((state[0] + 89.1 + INa_CaMh_half) / 6.086));
		double thsp = 3.0 * ths;
		state[13] = hssp - (hssp - state[13]) * exp(-dt / thsp);
		double hp = Ahf * state[10] + Ahs * state[13];
		double tjp = 1.46 * tj;
		state[14] = jss - (jss - state[14]) * exp(-dt / tjp);
		// printf("debug _ 1\n");
		double GNa = 75;
		double fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
		INa = GNa_ISO * GNa * (state[0] - ENa) * state[9] * state[9] * state[9] * ((1.0 - fINap) * h * state[12] + fINap * hp * state[14]);
	}


	/*TNNP INa formula*/

	/*double mss =1. / ((1. + exp((-56.86 - state[0]) / 9.03)) * (1. + exp((-56.86 - state[0]) / 9.03)));
	double tm = 1.0 / (6.765 * exp((state[0] + 11.64) / 34.77) + 8.552 * exp(-(state[0] + 77.42) / 5.955));
	state[9] = mss - (mss - state[9]) * exp(-dt / tm);


	double hss = 1. / ((1. + exp((state[0] + 71.55+ INah_half) / 7.43)) * (1. + exp((state[0] + 71.55+ INah_half) / 7.43)));
	double thf = 1.0 / (1.432e-5 * exp(-(state[0] + 1.196) / 6.285) + 6.149 * exp((state[0] + 0.5096) / 20.27));
	double ths = 1.0 / (0.009794 * exp(-(state[0] + 17.95) / 28.05) + 0.3343 * exp((state[0] + 5.730) / 56.66));
	double Ahf = 0.99;
	double Ahs = 1.0 - Ahf;
	state[10] = hss - (hss - state[10]) * exp(-dt / thf);
	state[11] = hss - (hss - state[11]) * exp(-dt / ths);
	double h = Ahf * state[10] + Ahs * state[11];
	double jss = hss;
	double tj = 2.038 + 1.0 / (0.02136 * exp(-(state[0] + 100.6) / 8.281) + 0.3052 * exp((state[0] + 0.9941) / 38.45));
	state[12] = jss - (jss - state[12]) * exp(-dt / tj);
	double hssp = 1.0 / (1 + exp((state[0] + 89.1 + INa_CaMh_half) / 6.086));
	double thsp = 3.0 * ths;
	state[13] = hssp - (hssp - state[13]) * exp(-dt / thsp);
	double hp = Ahf * state[10] + Ahs * state[13];
	double tjp = 1.46 * tj;
	state[14] = jss - (jss - state[14]) * exp(-dt / tjp);
	// printf("debug _ 1\n");
	double GNa = 75;
	double fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
	INa = GNa_ISO * GNa * (state[0] - ENa) * state[9] * state[9] * state[9] * ((1.0 - fINap) * h * state[12] + fINap * hp * state[14]);*/


	double tmL = tm;
	{	double mLss = 1.0 / (1.0 + exp((-(state[0] + 42.85)) / 5.264));
		state[15] = mLss - (mLss - state[15]) * exp(-dt / tmL);
		double hLss = 1.0 / (1.0 + exp((state[0] + 87.61) / 7.488));
		double thL = 200.0;
		state[16] = hLss - (hLss - state[16]) * exp(-dt / thL);
		double hLssp = 1.0 / (1.0 + exp((state[0] + 93.81) / 7.488));
		double thLp = 3.0 * thL;
		state[17] = hLssp - (hLssp - state[17]) * exp(-dt / thLp);
		double GNaL = 0.0075;
		if (celltype == LVEPI)
		{
			GNaL *= 0.6;
		}
		double fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
		INaL = GNaL * (state[0] - ENa) * state[15] * ((1.0 - fINaLp) * state[16] + fINaLp * state[17]);
	}

	{	double ass = 1.0 / (1.0 + exp((-(state[0] - 14.34)) / 14.82));
		double ta = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(state[0] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((state[0] + 100.0) / 29.3814)));
		state[18] = ass - (ass - state[18]) * exp(-dt / ta);
		double iss = 1.0 / (1.0 + exp((state[0] + 43.94 + InAcTo_Vhalf_ABh) / 5.711));
		double delta_LVEPI;
		if (celltype == LVEPI)
		{
			delta_LVEPI = 1.0 - (0.95 / (1.0 + exp((state[0] + 70.0) / 5.0)));
		}
		else
		{
			delta_LVEPI = 1.0;
		}
		double tiF = 4.562 + 1 / (0.3933 * exp((-(state[0] + 100.0)) / 100.0) + 0.08004 * exp((state[0] + 50.0) / 16.59));
		double tiS = 23.62 + 1 / (0.001416 * exp((-(state[0] + 96.52)) / 59.05) + 1.780e-8 * exp((state[0] + 114.1) / 8.079));
		tiF *= delta_LVEPI;
		tiS *= delta_LVEPI;
		double AiF = 1.0 / (1.0 + exp((state[0] - 213.6) / 151.2));
		double AiS = 1.0 - AiF;
		state[19] = iss - (iss - state[19]) * exp(-dt / tiF);
		state[20] = iss - (iss - state[20]) * exp(-dt / tiS);
		double i = AiF * state[19] + AiS * state[20];
		double assp = 1.0 / (1.0 + exp((-(state[0] - 24.34)) / 14.82));
		state[21] = assp - (assp - state[21]) * exp(-dt / ta);
		double dti_develop = 1.354 + 1.0e-4 / (exp((state[0] - 167.4) / 15.89) + exp(-(state[0] - 12.23) / 0.2154));
		double dti_recover = 1.0 - 0.5 / (1.0 + exp((state[0] + 70.0) / 20.0));
		double tiFp = dti_develop * dti_recover * tiF;
		double tiSp = dti_develop * dti_recover * tiS;
		state[22] = iss - (iss - state[22]) * exp(-dt / tiFp);
		state[23] = iss - (iss - state[23]) * exp(-dt / tiSp);
		double ip = AiF * state[22] + AiS * state[23];
		double Gto = 0.02;
		if (celltype == LVEPI)
		{
			Gto *= 4.0;
			// printf("LVEPI\n");
		}
		if (celltype == LVMCELL)
		{
			Gto *= 4.0;
			// printf("LVMCELL\n");
		}
		double fItop = (1.0 / (1.0 + KmCaMK / CaMKa));
		Ito = GTo_ABh * Gto * (state[0] - EK) * ((1.0 - fItop) * state[18] * i + fItop * state[21] * ip);
	}
	// printf("debug _ 2\n");

	{	double dss = 1.0 / (1.0 + exp((-(state[0] + 3.940 + ICaLd_half)) / 4.230));


		double td = 0.6 + 1.0 / (exp(-0.05 * (state[0] + 6.0)) + exp(0.09 * (state[0] + 14.0)));
		state[24] = dss - (dss - state[24]) * exp(-dt / td);

		// double fss = 1.0 / (1.0 + exp((state[0] + 19.58 + ICaLf_half) / 3.696));  // ref 1
		// double fss = 1.0 / (1.0 + exp((state[0] + 20.5 + ICaLf_half) / 3.1775));  // ref 1
		double fss = 1.0 / (1.0 + exp((state[0] + 20.25 + ICaLf_half) / 3.696));  // ref 1

		// double fss = 1.0 / (1.0 + exp((state[0] + 19.58 + ICaLf_half) / 3.1775));  // ref 1
		// double fss = (1.0 / (1.0 + exp((state[0] + 21.979 + ICaLf_half) / 3.1775)) + 0.0001) * (1.0 / 1.0001); // ref 2
		// double fss = (1.0 / (1.0 + exp((state[0] + 20.979 + ICaLf_half) / 3.1775)) ) ; // ref 2




		double tff = Tau_ICaL_ff_ISO * (7.0 + 1.0 / (0.0045 * exp(-(state[0] + 20.0) / 10.0) + 0.0045 * exp((state[0] + 20.0) / 10.0)));
		double tfs = 1000.0 + 1.0 / (0.000035 * exp(-(state[0] + 5.0) / 4.0) + 0.000035 * exp((state[0] + 5.0) / 6.0));
		double Aff = 0.6;
		double Afs = 1.0 - Aff;
		state[25] = fss - (fss - state[25]) * exp(-dt / tff);
		state[26] = fss - (fss - state[26]) * exp(-dt / tfs);
		double f = Aff * state[25] + Afs * state[26];
		double fcass = fss;
		double tfcaf = 7.0 + 1.0 / (0.04 * exp(-(state[0] - 4.0) / 7.0) + 0.04 * exp((state[0] - 4.0) / 7.0));
		double tfcas = 100.0 + 1.0 / (0.00012 * exp(-state[0] / 3.0) + 0.00012 * exp(state[0] / 7.0));
		double Afcaf = 0.3 + 0.6 / (1.0 + exp((state[0] - 10.0) / 10.0));
		double Afcas = 1.0 - Afcaf;
		state[27] = fcass - (fcass - state[27]) * exp(-dt / tfcaf);
		state[28] = fcass - (fcass - state[28]) * exp(-dt / tfcas);
		double fca = Afcaf * state[27] + Afcas * state[28];
		double tjca = 75.0;
		state[29] = fcass - (fcass - state[29]) * exp(-dt / tjca);
		double tffp = 2.5 * tff;
		state[31] = fss - (fss - state[31]) * exp(-dt / tffp);
		double fp = Aff * state[31] + Afs * state[26];
		double tfcafp = 2.5 * tfcaf;
		state[32] = fcass - (fcass - state[32]) * exp(-dt / tfcafp);
		double fcap = Afcaf * state[32] + Afcas * state[28];
		double Kmn = 0.002;
		double k2n = 1000.0;
		double km2n = state[29] * 1.0;
		double anca = 1.0 / (k2n / km2n + pow(1.0 + Kmn / state[6], 4.0));
		state[30] = anca * k2n / km2n - (anca * k2n / km2n - state[30]) * exp(-km2n * dt);

		double Cass = state[6];

		if (Cass > 0.03)
		{
			Cass = 0.03;   // 0.03mM ceiling, from ref 1. * 1. Ohara et al. Arrhythmia formation in subclinical (“silent”) long QT syndrome requires multiple insults: Quantitative mechanistic study using the KCNQ1 mutation Q357R as example
		}
		double PhiCaL = 4.0 * vffrt * (Cass * exp(2.0 * vfrt) - 0.341 * cao) / (exp(2.0 * vfrt) - 1.0);
		double PhiCaNa = 1.0 * vffrt * (0.75 * state[2] * exp(1.0 * vfrt) - 0.75 * nao) / (exp(1.0 * vfrt) - 1.0);
		double PhiCaK = 1.0 * vffrt * (0.75 * state[4] * exp(1.0 * vfrt) - 0.75 * ko) / (exp(1.0 * vfrt) - 1.0);

		double PCa = Pca_ISO * 0.0001;
		if (celltype == LVEPI)
		{
			PCa *= 1.2;
		}
		if (celltype == LVMCELL)
		{
			PCa *= 2.5;
		}
		double PCap = 1.1 * PCa;
		double PCaNa = 0.00125 * PCa;
		double PCaK = 3.574e-4 * PCa;
		double PCaNap = 0.00125 * PCap;
		double PCaKp = 3.574e-4 * PCap;
		double fICaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
		ICaL = (1.0 - fICaLp) * PCa * PhiCaL * state[24] * (f * (1.0 - state[30]) + state[29] * fca * state[30]) + fICaLp * PCap * PhiCaL * state[24] * (fp * (1.0 - state[30]) + state[29] * fcap * state[30]);
		ICaNa = (1.0 - fICaLp) * PCaNa * PhiCaNa * state[24] * (f * (1.0 - state[30]) + state[29] * fca * state[30]) + fICaLp * PCaNap * PhiCaNa * state[24] * (fp * (1.0 - state[30]) + state[29] * fcap * state[30]);
		ICaK = (1.0 - fICaLp) * PCaK * PhiCaK * state[24] * (f * (1.0 - state[30]) + state[29] * fca * state[30]) + fICaLp * PCaKp * PhiCaK * state[24] * (fp * (1.0 - state[30]) + state[29] * fcap * state[30]);
	}
	{	double xrss = 1.0 / (1.0 + exp((-(state[0] + 8.337)) / 6.789));
		double txrf = 12.98 + 1.0 / (0.3652 * exp((state[0] - 31.66) / 3.869) + 4.123e-5 * exp((-(state[0] - 47.78)) / 20.38));
		double txrs = 1.865 + 1.0 / (0.06629 * exp((state[0] - 34.70) / 7.355) + 1.128e-5 * exp((-(state[0] - 29.74)) / 25.94));
		double Axrf = 1.0 / (1.0 + exp((state[0] + 54.81) / 38.21));
		double Axrs = 1.0 - Axrf;
		state[33] = xrss - (xrss - state[33]) * exp(-dt / txrf);
		state[34] = xrss - (xrss - state[34]) * exp(-dt / txrs);
		double xr = Axrf * state[33] + Axrs * state[34];
		double rkr = 1.0 / (1.0 + exp((state[0] + 55.0) / 75.0)) * 1.0 / (1.0 + exp((state[0] - 10.0) / 30.0));
		double GKr = 0.046;
		if (celltype == LVEPI)
		{
			GKr *= 1.3;
		}
		if (celltype == LVMCELL)
		{
			// GKr *= 0.8;
			GKr *= 0.8;  // 1.3*0.91 in later paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3311457/
		}
		IKr = GKr * sqrt(ko / 5.4) * xr * rkr * (state[0] - EK);
	}


	if (mutation == ORd)
	{
		double xs1ss = 1.0 / (1.0 + exp((-(state[0] + IKs_ac_shift_ISO + 11.60)) / 8.932));
		double txs1 = Tau_Ks * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;
		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.0193 * exp((-(state[0] + 66.54)) / 31.0));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034;
		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);
	} else if (mutation == IKS_WT)
	{

		double xs1ss = 1.0 / (1.0 + exp((-(state[0] + 2.4 + IKs_ac_shift_ISO)) / 13.2));
		// in the new IKs model, txs1 was timed by 0.5 to better reproduce the results of voltage clamp.

		double txs1 = Tau_Ks * 0.5 * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;


		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.0071038 * exp((-(state[0] + 90.5741763432)) / 69.080));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034 * 1.2; // to make sure that the new IKS_WT IKs produced same APD90 as the original ORd IKs model
		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);
	} else if (mutation == HET)
	{
		double xs1ss = 1.0 / (1.0 + exp((-(state[0] - 25.4 + IKs_ac_shift_ISO)) / 17.8));
		double txs1 = Tau_Ks * 0.5 * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;


		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.017671 * exp((-(state[0] + 93.997 )) / 43.8948649349));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034 * 1.2; // to make sure that the new IKS_WT IKs produced same APD90 as the original ORd IKs model,
		// and the conductance of IKs has been scaled so that the tail current at 30mV in the voltage clamp matches ratio HET/IKS_WT in the original experimetnal data
		// see Figure 2E in paper Wu et al. 2014. A Molecular Mechanism for Adrenergic-Induced Long QT Syndrome  Jie Wu, P H D,*yz Nobu Naiki, MD,y Wei-Guang Ding, MD, P H D,z

		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);

	} else if (mutation == G269S)
	{
		double xs1ss = 1.0 / (1.0 + exp((-(state[0] - 68.3 + IKs_ac_shift_ISO)) / 23.4));
		double txs1 = Tau_Ks * 0.5 * ( 817.3 + 1.0 / (2.326e-4 * exp((state[0] + 48.28) / 17.80) + 0.001292 * exp((-(state[0] + 210.0)) / 230.0)));
		state[35] = xs1ss - (xs1ss - state[35]) * exp(-dt / txs1);
		double xs2ss = xs1ss;

		double txs2 = 1.0 / (0.01 * exp((state[0] - 50.0) / 20.0) + 0.0290693 * exp((-(state[0] +  106.263263059)) / 48.7397674066));
		state[36] = xs2ss - (xs2ss - state[36]) * exp(-dt / txs2);
		double KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-5 / state[5], 1.4));
		double GKs = 0.0034 * 5.2 * 1.2;
		if (celltype == LVEPI)
		{
			GKs *= 1.4;
		}
		IKs = GKs_ABh * GKs_ISO *  GKs * KsCa * state[35] * state[36] * (state[0] - EKs);
	} else {
		printf("wrong IKs mutation option!!!\n");
		exit(0);
	}

	double xk1ss = 1.0 / (1.0 + exp(-(state[0] + 2.5538 * ko + 144.59) / (1.5692 * ko + 3.8115)));
	double txk1 = 122.2 / (exp((-(state[0] + 127.2)) / 20.36) + exp((state[0] + 236.8) / 69.33));
	state[37] = xk1ss - (xk1ss - state[37]) * exp(-dt / txk1);
	double rk1 = 1.0 / (1.0 + exp((state[0] + 105.8 - 2.6 * ko) / 9.493));
	double GK1 = 0.1908;
	if (celltype == LVEPI)
	{
		GK1 *= 1.2;
	}
	if (celltype == LVMCELL)
	{
		GK1 *= 1.3;
	}
	IK1 = GK1 * sqrt(ko) * rk1 * state[37] * (state[0] - EK);

	// double 15.0 = 15.0;
	// double 5.0 = 5.0;
	// double 88.12 = 88.12;
	// double 12.5 = 12.5;
	double x1;
	double x2;
	double x3;
	double x4;
	double E1, E2, E3, E4;
	{	double wna = 6.0e4;
		double wca = 6.0e4;
		double wnaca = 5.0e3;
		double kcaon = 1.5e6;
		double kcaoff = 5.0e3;
		double qna = 0.5224;
		double qca = 0.1670;
		double hca = exp((qca * state[0] * F) / (R * T));
		double hna = exp((qna * state[0] * F) / (R * T));
		double h1 = 1 + state[1] / 88.12 * (1 + hna);
		double h2 = (state[1] * hna) / (88.12 * h1);
		double h3 = 1.0 / h1;
		double h4 = 1.0 + state[1] / 15.0 * (1 + state[1] / 5.0);
		double h5 = state[1] * state[1] / (h4 * 15.0 * 5.0);
		double h6 = 1.0 / h4;
		double h7 = 1.0 + nao / 88.12 * (1.0 + 1.0 / hna);
		double h8 = nao / (88.12 * hna * h7);
		double h9 = 1.0 / h7;
		double h10 = 12.5 + 1.0 + nao / 15.0 * (1.0 + nao / 5.0);
		double h11 = nao * nao / (h10 * 15.0 * 5.0);
		double h12 = 1.0 / h10;
		double k1 = h12 * cao * kcaon;
		double k2 = kcaoff;
		double k3p = h9 * wca;
		double k3pp = h8 * wnaca;
		double k3 = k3p + k3pp;
		double k4p = h3 * wca / hca;
		double k4pp = h2 * wnaca;
		double k4 = k4p + k4pp;
		double k5 = kcaoff;
		double k6 = h6 * state[5] * kcaon;
		double k7 = h5 * h2 * wna;
		double k8 = h8 * h11 * wna;
		x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
		x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
		x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
		x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
		E1 = x1 / (x1 + x2 + x3 + x4);
		E2 = x2 / (x1 + x2 + x3 + x4);
		E3 = x3 / (x1 + x2 + x3 + x4);
		E4 = x4 / (x1 + x2 + x3 + x4);
		double KmCaAct = 150.0e-6;
		double allo = 1.0 / (1.0 + pow(KmCaAct / state[5], 2.0));

		double JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
		double JncxCa = E2 * k2 - E1 * k1;
		double Gncx = 0.0008;
		if (celltype == LVEPI)
		{
			Gncx *= 1.1;
		}
		if (celltype == LVMCELL)
		{
			Gncx *= 1.4;
		}
		INaCa_i = 0.8 * Gncx * allo * (zna * JncxNa + zca * JncxCa);

		h1 = 1 + state[2] / 88.12 * (1 + hna);
		h2 = (state[2] * hna) / (88.12 * h1);
		h3 = 1.0 / h1;
		h4 = 1.0 + state[2] / 15.0 * (1 + state[2] / 5.0);
		h5 = state[2] * state[2] / (h4 * 15.0 * 5.0);
		h6 = 1.0 / h4;
		h7 = 1.0 + nao / 88.12 * (1.0 + 1.0 / hna);
		h8 = nao / (88.12 * hna * h7);
		h9 = 1.0 / h7;
		h10 = 12.5 + 1.0 + nao / 15.0 * (1 + nao / 5.0);
		h11 = nao * nao / (h10 * 15.0 * 5.0);
		h12 = 1.0 / h10;
		k1 = h12 * cao * kcaon;
		k2 = kcaoff;
		k3p = h9 * wca;
		k3pp = h8 * wnaca;
		k3 = k3p + k3pp;
		k4p = h3 * wca / hca;
		k4pp = h2 * wnaca;
		k4 = k4p + k4pp;
		k5 = kcaoff;
		k6 = h6 * state[6] * kcaon;
		k7 = h5 * h2 * wna;
		k8 = h8 * h11 * wna;
		x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
		x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
		x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
		x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
		E1 = x1 / (x1 + x2 + x3 + x4);
		E2 = x2 / (x1 + x2 + x3 + x4);
		E3 = x3 / (x1 + x2 + x3 + x4);
		E4 = x4 / (x1 + x2 + x3 + x4);
		KmCaAct = 150.0e-6;
		allo = 1.0 / (1.0 + pow(KmCaAct / state[6], 2.0));
		JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
		JncxCa = E2 * k2 - E1 * k1;
		INaCa_ss = 0.2 * Gncx * allo * (zna * JncxNa + zca * JncxCa);

		INaCa = INaCa_i + INaCa_ss;
	}
	// printf("debug _ 3\n");

	{	double k1p = 949.5;
		double k1m = 182.4;
		double k2p = 687.2;
		double k2m = 39.4;
		double k3p = 1899.0;
		double k3m = 79300.0;
		double k4p = 639.0;
		double k4m = 40.0;
		// double Knai0 = 9.073;
		// double Knao0 = 27.78;
		double delta = -0.1550;
		double Knai = NaK_ISO * (9.073 * exp((delta * state[0] * F) / (3.0 * R * T)));
		double Knao = 27.78 * exp(((1.0 - delta) * state[0] * F) / (3.0 * R * T));
		double Kki = 0.5;
		double Kko = 0.3582;
		double MgADP = 0.05;
		double MgATP = 9.8;
		double Kmgatp = 1.698e-7;
		double H = 1.0e-7;
		double eP = 4.2;
		double Khp = 1.698e-7;
		double Knap = 224.0;
		double Kxkur = 292.0;
		double P = eP / (1.0 + H / Khp + state[1] / Knap + state[3] / Kxkur);
		double a1 = (k1p * pow(state[1] / Knai, 3.0)) / (pow(1.0 + state[1] / Knai, 3.0) + pow(1.0 + state[3] / Kki, 2.0) - 1.0);
		double b1 = k1m * MgADP;
		double a2 = k2p;
		double b2 = (k2m * pow(nao / Knao, 3.0)) / (pow(1.0 + nao / Knao, 3.0) + pow(1.0 + ko / Kko, 2.0) - 1.0);
		double a3 = (k3p * pow(ko / Kko, 2.0)) / (pow(1.0 + nao / Knao, 3.0) + pow(1.0 + ko / Kko, 2.0) - 1.0);
		double b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp);
		double a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp);
		double b4 = (k4m * pow(state[3] / Kki, 2.0)) / (pow(1.0 + state[1] / Knai, 3.0) + pow(1.0 + state[3] / Kki, 2.0) - 1.0);
		double x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
		double x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
		double x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
		double x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
		E1 = x1 / (x1 + x2 + x3 + x4);
		E2 = x2 / (x1 + x2 + x3 + x4);
		E3 = x3 / (x1 + x2 + x3 + x4);
		E4 = x4 / (x1 + x2 + x3 + x4);
		double zk = 1.0;
		double JnakNa = 3.0 * (E1 * a3 - E2 * b3);
		double JnakK = 2.0 * (E4 * b1 - E3 * a1);
		double Pnak = 30;
		if (celltype == LVEPI)
		{
			Pnak *= 0.9;
		}
		if (celltype == LVMCELL)
		{
			Pnak *= 0.7;
		}
		INaK = Pnak * (zna * JnakNa + zk * JnakK);

		double xkb = 1.0 / (1.0 + exp(-(state[0] - 14.48) / 18.34));
		double GKb = 0.003;
		IKb = GKb_ISO * GKb * xkb * (state[0] - EK);

		double PNab = 3.75e-10;
		INab = PNab * vffrt * (state[1] * exp(vfrt) - nao) / (exp(vfrt) - 1.0);

		double PCab = 2.5e-8;
		ICab = PCab * 4.0 * vffrt * (state[5] * exp(2.0 * vfrt) - 0.341 * cao) / (exp(2.0 * vfrt) - 1.0);

		/*double GpCa = 0.005;
		IpCa = GpCa * state[5] / (0.00025 + state[5]);  // note the cai pump was increased to expel more cai.*/
		double GpCa = 0.0005;
		IpCa = GpCa * state[5] / (0.0005 + state[5]);
	}


	// // original..
	//  double GpCa = 0.0005;
	// IpCa = GpCa * cai / (0.0005 + cai);

	JdiffNa = (state[2] - state[1]) / 2.0;
	JdiffK = (state[4] - state[3]) / 2.0;
	Jdiff = (state[6] - state[5]) / 0.2;

	double bt = 4.75;
	double a_rel = A_rel_ISO  * 0.5 * bt ;

	double Jrel_inf = a_rel * (-ICaL) / (1.0 + pow(1.5 / state[8], 8.0));
	if (celltype == LVMCELL)
	{
		Jrel_inf *= 1.7;
	}
	double tau_rel = tau_rel_ISO * bt / (1.0 + 0.0123 / state[8]);
	if (tau_rel < 0.005)
	{
		tau_rel = 0.005;
	}
	state[38] = Jrel_inf - (Jrel_inf - state[38]) * exp(-dt / tau_rel);
	double btp = 1.25 * bt;
	double a_relp = A_rel_ISO * 0.5 * btp;
	double Jrel_infp = a_relp * (-ICaL) / (1.0 + pow(1.5 / state[8], 8.0));
	if (celltype == LVMCELL)
	{
		Jrel_infp *= 1.7;
	}
	double tau_relp = tau_rel_ISO * btp / (1.0 + 0.0123 / state[8]);
	if (tau_relp < 0.005)
	{
		tau_relp = 0.005;
	}
	state[39] = Jrel_infp - (Jrel_infp - state[39]) * exp(-dt / tau_relp);
	double fJrelp = (1.0 / (1.0 + KmCaMK / CaMKa));
	Jrel = (1.0 - fJrelp) * state[38] + fJrelp * state[39];




	double Jupnp = 0.004375 * state[5] / (state[5] + 0.00092 * Jup_ISO);
	double Jupp = 2.75 * 0.004375 * state[5] / (state[5] + ( 0.00092 - 0.00017) * Jup_ISO);


	if (celltype == LVEPI)
	{
		Jupnp *= 1.3;
		Jupp *= 1.3;
	}


	double fJupp = (1.0 / (1.0 + KmCaMK / CaMKa));
	Jleak = 0.0039375 * state[7] / 15.0;
	Jup = (1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak;

	Jtr = (state[7] - state[8]) / 100.0;

	state[1] += dt * (-(INa + INaL + 3.0 * INaCa_i + 3.0 * INaK + INab) * Acap / (F * vmyo) + JdiffNa * vss / vmyo);
	state[2] += dt * (-(ICaNa + 3.0 * INaCa_ss) * Acap / (F * vss) - JdiffNa);

	state[3] += dt * (-(Ito + IKr + IKs + IK1 + IKb + Ist - 2.0 * INaK) * Acap / (F * vmyo) + JdiffK * vss / vmyo);
	state[4] += dt * (-(ICaK) * Acap / (F * vss) - JdiffK);
	double kmtrpn_cell = kmtrpn * kmtrpn_ISO;
	double Bcai;
	if (celltype == LVEPI)
	{
		Bcai = 1.0 / (1.0 + 1.3 * cmdnmax * kmcmdn / pow(kmcmdn + state[5], 2.0) + trpnmax * kmtrpn_cell / pow(kmtrpn_cell + state[5], 2.0));
	}
	else
	{
		Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / pow(kmcmdn + state[5], 2.0) + trpnmax * kmtrpn_cell / pow(kmtrpn_cell + state[5], 2.0));
	}
	state[5] += dt * (Bcai * (-(IpCa + ICab - 2.0 * INaCa_i) * Acap / (2.0 * F * vmyo) - Jup * vnsr / vmyo + Jdiff * vss / vmyo));

	double Bcass = 1.0 / (1.0 + BSRmax * KmBSR / pow(KmBSR + state[6], 2.0) + BSLmax * KmBSL / pow(KmBSL + state[6], 2.0));
	state[6] += dt * (Bcass * (-(ICaL - 2.0 * INaCa_ss) * Acap / (2.0 * F * vss) + Jrel * vjsr / vss - Jdiff));

	state[7] += dt * (Jup - Jtr * vjsr / vnsr);

	double Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / pow(kmcsqn + state[8], 2.0));
	state[8] += dt * (Bcajsr * (Jtr - Jrel));

	return  (INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + IKb + IpCa + ICab + Ist);
	// f// printf(stderr, "Ina2=%f\n", INa);

	// double I = INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + IKb + IpCa + ICab + Ist;
	// f// printf(stderr, "I=%f\n", I);
	// printf("debug _ 4\n");

	// ion -> IKr = IKr;
	// ion ->  IKs = IKs;
	// ion -> IK1 = IK1;
	// ion -> Ito = Ito;
	// // ion -> Isus = Isus;
	// ion -> INa = INa;
	// ion ->  IbNa = INab;
	// ion ->ICaL =  ICaL;
	// ion -> IbCa = ICab;
	// ion -> INaK = INaK;
	// ion ->  INaCa = INaCa;//!!!!!
	// ion -> IpCa = IpCa;
	// ion -> Cai = state[5];//!!!!

	// ion -> IpK = IpK;
	// ion -> If = If;
	// printf("debug _ 5\n");


	/*static double t;
	t += dt;
	static int count;
	count += 1;
	FILE * output = fopen("current.txt", "a");
	if (count % 10 == 0)
	{

	    fprintf(output, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t,Jup, Jrel, INa , INaL , Ito , ICaL , ICaNa , ICaK , IKr , IKs , IK1 , INaCa , INaK , INab , IKb , IpCa , ICab , Ist);
	}
	fclose(output);*/

	// IKs_rec = IKs;


}

void initial_rudy_state(double *state, TypeCell celtype)
{
	state[0] = -87.5;
	state[1] = 7;
	state[2] = 7;
	state[3] = 145;
	state[4] = 145;
	state[5] = 1.0e-4;
	state[6] = 1.0e-4;
	state[7] = 1.2;
	state[8] = 1.2;
	state[9] = 0;
	state[10] = 1;
	state[11] = 1;
	state[12] = 1;
	state[13] = 1;
	state[14] = 1;
	state[15] = 0;
	state[16] = 1;
	state[17] = 1;
	state[18] = 0;
	state[19] = 1;
	state[20] = 1;
	state[21] = 0;
	state[22] = 1;
	state[23] = 1;
	state[24] = 0;
	state[25] = 1;
	state[26] = 1;
	state[27] = 1;
	state[28] = 1;
	state[29] = 1;
	state[30] = 0;
	state[31] = 1;
	state[32] = 1;
	state[33] = 0;
	state[34] = 0;
	state[35] = 0;
	state[36] = 0;
	state[37] = 1;
	state[38] = 0;
	state[39] = 0;
	state[40] = 0;
}


void INa_Drug_initial_rudy_state(double *state, TypeCell celtype)
{
	state[0] = -87.5;
	state[1] = 7;
	state[2] = 7;
	state[3] = 145;
	state[4] = 145;
	state[5] = 1.0e-4;
	state[6] = 1.0e-4;
	state[7] = 1.2;
	state[8] = 1.2;
	state[9] = 0;
	state[10] = 1;
	state[11] = 1;
	state[12] = 1;
	state[13] = 1;
	state[14] = 1;
	state[15] = 0;
	state[16] = 1;
	state[17] = 1;
	state[18] = 0;
	state[19] = 1;
	state[20] = 1;
	state[21] = 0;
	state[22] = 1;
	state[23] = 1;
	state[24] = 0;
	state[25] = 1;
	state[26] = 1;
	state[27] = 1;
	state[28] = 1;
	state[29] = 1;
	state[30] = 0;
	state[31] = 1;
	state[32] = 1;
	state[33] = 0;
	state[34] = 0;
	state[35] = 0;
	state[36] = 0;
	state[37] = 1;
	state[38] = 0;
	state[39] = 0;
	state[40] = 0;
	state[41] = 0;
	state[42] = 0;

}



int ORd_con_length() {
	return 43;
}