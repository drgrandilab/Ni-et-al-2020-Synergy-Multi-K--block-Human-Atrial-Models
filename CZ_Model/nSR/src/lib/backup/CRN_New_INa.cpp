#include <cmath>
#include "CRN_ODE.h"
#include "ExplicitSolver.hpp"


double CRN_1998_ODE_New_INa(const double dt, const double Istim,  SingleCellPara & data, double *state) {

	//------------------------------------------------------------------------------
	// Constants
	//------------------------------------------------------------------------------
	double CMDN_max;   // millimolar (in Ca_buffers)
	double CSQN_max;   // millimolar (in Ca_buffers)
	double Km_CMDN;   // millimolar (in Ca_buffers)
	double Km_CSQN;   // millimolar (in Ca_buffers)
	double Km_TRPN;   // millimolar (in Ca_buffers)
	double TRPN_max;   // millimolar (in Ca_buffers)
	double Ca_up_max;   // millimolar (in Ca_leak_current_by_the_NSR)
	double K_rel;   // per_millisecond (in Ca_release_current_from_JSR)
	double I_up_max;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
	double K_up;   // millimolar (in Ca_uptake_current_by_the_NSR)
	double g_Ca_L;   // nanoS_per_picoF (in L_type_Ca_channel)
	double I_NaCa_max;   // picoA_per_picoF (in Na_Ca_exchanger_current)
	double K_mCa;   // millimolar (in Na_Ca_exchanger_current)
	double K_mNa;   // millimolar (in Na_Ca_exchanger_current)
	double K_sat;   // dimensionless (in Na_Ca_exchanger_current)
	double gamma;   // dimensionless (in Na_Ca_exchanger_current)
	double g_B_Ca;   // nanoS_per_picoF (in background_currents)
	double g_B_K;   // nanoS_per_picoF (in background_currents)
	double g_B_Na;   // nanoS_per_picoF (in background_currents)
	double g_Na;   // nanoS_per_picoF (in fast_sodium_current)
	double V_cell;   // micrometre_3 (in intracellular_ion_concentrations)
	double Cm;   // picoF (in membrane)
	double F;   // coulomb_per_millimole (in membrane)
	double R;   // joule_per_mole_kelvin (in membrane)
	double T;   // kelvin (in membrane)
	double stim_amplitude;   // picoA (in membrane)
	double stim_duration;   // millisecond (in membrane)
	double stim_end;   // millisecond (in membrane)
	double stim_period;   // millisecond (in membrane)
	double stim_start;   // millisecond (in membrane)
	double g_Kr;   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
	double i_CaP_max;   // picoA_per_picoF (in sarcolemmal_calcium_pump_current)
	double g_Ks;   // nanoS_per_picoF (in slow_delayed_rectifier_K_current)
	double Km_K_o;   // millimolar (in sodium_potassium_pump)
	double Km_Na_i;   // millimolar (in sodium_potassium_pump)
	double i_NaK_max;   // picoA_per_picoF (in sodium_potassium_pump)
	double Ca_o;   // millimolar (in standard_ionic_concentrations)
	double K_o;   // millimolar (in standard_ionic_concentrations)
	double Na_o;   // millimolar (in standard_ionic_concentrations)
	double g_K1;   // nanoS_per_picoF (in time_independent_potassium_current)
	double tau_tr;   // millisecond (in transfer_current_from_NSR_to_JSR)
	double K_Q10;   // dimensionless (in transient_outward_K_current)
	double g_to;   // nanoS_per_picoF (in transient_outward_K_current)

	//------------------------------------------------------------------------------
	// Computed variables
	//------------------------------------------------------------------------------

	double Ca_CMDN;   // millimolar (in Ca_buffers)
	double Ca_CSQN;   // millimolar (in Ca_buffers)
	double Ca_TRPN;   // millimolar (in Ca_buffers)
	double i_up_leak;   // millimolar_per_millisecond (in Ca_leak_current_by_the_NSR)
	double tau_u;   // millisecond (in Ca_release_current_from_JSR_u_gate)
	double u_infinity;   // dimensionless (in Ca_release_current_from_JSR_u_gate)
	double tau_v;   // millisecond (in Ca_release_current_from_JSR_v_gate)
	double v_infinity;   // dimensionless (in Ca_release_current_from_JSR_v_gate)
	double tau_w;   // millisecond (in Ca_release_current_from_JSR_w_gate)
	double w_infinity;   // dimensionless (in Ca_release_current_from_JSR_w_gate)
	double Fn;   // dimensionless (in Ca_release_current_from_JSR)
	double i_rel;   // millimolar_per_millisecond (in Ca_release_current_from_JSR)
	double i_up;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
	double d_infinity;   // dimensionless (in L_type_Ca_channel_d_gate)
	double tau_d;   // millisecond (in L_type_Ca_channel_d_gate)
	double f_Ca_infinity;   // dimensionless (in L_type_Ca_channel_f_Ca_gate)
	double tau_f_Ca;   // millisecond (in L_type_Ca_channel_f_Ca_gate)
	double f_infinity;   // dimensionless (in L_type_Ca_channel_f_gate)
	double tau_f;   // millisecond (in L_type_Ca_channel_f_gate)
	double i_Ca_L;   // picoA (in L_type_Ca_channel)
	double i_NaCa;   // picoA (in Na_Ca_exchanger_current)
	double E_Ca;   // millivolt (in background_currents)
	double i_B_Ca;   // picoA (in background_currents)
	double i_B_K;   // picoA (in background_currents)
	double i_B_Na;   // picoA (in background_currents)
	double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
	double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
	double h_inf;   // dimensionless (in fast_sodium_current_h_gate)
	double tau_h;   // millisecond (in fast_sodium_current_h_gate)
	double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
	double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
	double j_inf;   // dimensionless (in fast_sodium_current_j_gate)
	double tau_j;   // millisecond (in fast_sodium_current_j_gate)
	double alpha_m;   // per_millisecond (in fast_sodium_current_m_gate)
	double beta_m;   // per_millisecond (in fast_sodium_current_m_gate)
	double m_inf;   // dimensionless (in fast_sodium_current_m_gate)
	double tau_m;   // millisecond (in fast_sodium_current_m_gate)
	double E_Na;   // millivolt (in fast_sodium_current)
	double i_Na;   // picoA (in fast_sodium_current)
	double B1;   // millimolar_per_millisecond (in intracellular_ion_concentrations)
	double B2;   // dimensionless (in intracellular_ion_concentrations)
	double V_i;   // micrometre_3 (in intracellular_ion_concentrations)
	double V_rel;   // micrometre_3 (in intracellular_ion_concentrations)
	double V_up;   // micrometre_3 (in intracellular_ion_concentrations)
	double i_st;   // picoA (in membrane)
	double alpha_xr;   // per_millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
	double beta_xr;   // per_millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
	double tau_xr;   // millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
	double xr_infinity;   // dimensionless (in rapid_delayed_rectifier_K_current_xr_gate)
	double i_Kr;   // picoA (in rapid_delayed_rectifier_K_current)
	double i_CaP;   // picoA (in sarcolemmal_calcium_pump_current)
	double alpha_xs;   // per_millisecond (in slow_delayed_rectifier_K_current_xs_gate)
	double beta_xs;   // per_millisecond (in slow_delayed_rectifier_K_current_xs_gate)
	double tau_xs;   // millisecond (in slow_delayed_rectifier_K_current_xs_gate)
	double xs_infinity;   // dimensionless (in slow_delayed_rectifier_K_current_xs_gate)
	double i_Ks;   // picoA (in slow_delayed_rectifier_K_current)
	double f_NaK;   // dimensionless (in sodium_potassium_pump)
	double i_NaK;   // picoA (in sodium_potassium_pump)
	double sigma;   // dimensionless (in sodium_potassium_pump)
	double E_K;   // millivolt (in time_independent_potassium_current)
	double i_K1;   // picoA (in time_independent_potassium_current)
	double i_tr;   // millimolar_per_millisecond (in transfer_current_from_NSR_to_JSR)
	double alpha_oa;   // per_millisecond (in transient_outward_K_current_oa_gate)
	double beta_oa;   // per_millisecond (in transient_outward_K_current_oa_gate)
	double oa_infinity;   // dimensionless (in transient_outward_K_current_oa_gate)
	double tau_oa;   // millisecond (in transient_outward_K_current_oa_gate)
	double alpha_oi;   // per_millisecond (in transient_outward_K_current_oi_gate)
	double beta_oi;   // per_millisecond (in transient_outward_K_current_oi_gate)
	double oi_infinity;   // dimensionless (in transient_outward_K_current_oi_gate)
	double tau_oi;   // millisecond (in transient_outward_K_current_oi_gate)
	double i_to;   // picoA (in transient_outward_K_current)
	double alpha_ua;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double beta_ua;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double tau_ua;   // millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double ua_infinity;   // dimensionless (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double alpha_ui;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double beta_ui;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double tau_ui;   // millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double ui_infinity;   // dimensionless (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double g_Kur;   // nanoS_per_picoF (in ultrarapid_delayed_rectifier_K_current)
	double i_Kur;   // picoA (in ultrarapid_delayed_rectifier_K_current)
	double phi;


	//---------------------------------------------------------------------------
	// Constants
	//---------------------------------------------------------------------------

	CMDN_max = 0.05;   // millimolar (in Ca_buffers)
	CSQN_max = 10.0;   // millimolar (in Ca_buffers)
	Km_CMDN = 0.00238;   // millimolar (in Ca_buffers)
	Km_CSQN = 0.8;   // millimolar (in Ca_buffers)
	Km_TRPN = 0.0005;   // millimolar (in Ca_buffers)
	TRPN_max = 0.07;   // millimolar (in Ca_buffers)
	Ca_up_max = 15.0;   // millimolar (in Ca_leak_current_by_the_NSR)
	K_rel = 30.0;   // per_millisecond (in Ca_release_current_from_JSR)
	I_up_max = 0.005;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
	K_up = 0.00092;   // millimolar (in Ca_uptake_current_by_the_NSR)
	g_Ca_L = 0.12375;   // nanoS_per_picoF (in L_type_Ca_channel)
	I_NaCa_max = 1600.0;   // picoA_per_picoF (in Na_Ca_exchanger_current)
	K_mCa = 1.38;   // millimolar (in Na_Ca_exchanger_current)
	K_mNa = 87.5;   // millimolar (in Na_Ca_exchanger_current)
	K_sat = 0.1;   // dimensionless (in Na_Ca_exchanger_current)
	gamma = 0.35;   // dimensionless (in Na_Ca_exchanger_current)
	g_B_Ca = 0.001131;   // nanoS_per_picoF (in background_currents)
	g_B_K = 0.0;   // nanoS_per_picoF (in background_currents)
	g_B_Na = 0.0006744375;   // nanoS_per_picoF (in background_currents)
	g_Na = 7.8;   // nanoS_per_picoF (in fast_sodium_current)
	V_cell = 20100.0;   // micrometre_3 (in intracellular_ion_concentrations)
	Cm = 100.0;   // picoF (in membrane)
	F = 96.4867;   // coulomb_per_millimole (in membrane)
	R = 8.3143;   // joule_per_mole_kelvin (in membrane)
	T = 310.0;   // kelvin (in membrane)
	stim_amplitude = -2000.0;   // picoA (in membrane)
	stim_duration = 2.0;   // millisecond (in membrane)
	stim_end = 50000.0;   // millisecond (in membrane)
	stim_period = 1000.0;   // millisecond (in membrane)
	stim_start = 50.0;   // millisecond (in membrane)
	g_Kr = 0.029411765;   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
	i_CaP_max = 0.275;   // picoA_per_picoF (in sarcolemmal_calcium_pump_current)
	g_Ks = 0.12941176;   // nanoS_per_picoF (in slow_delayed_rectifier_K_current)
	Km_K_o = 1.5;   // millimolar (in sodium_potassium_pump)
	Km_Na_i = 10.0;   // millimolar (in sodium_potassium_pump)
	i_NaK_max = 0.59933874;   // picoA_per_picoF (in sodium_potassium_pump)
	Ca_o = 1.8;   // millimolar (in standard_ionic_concentrations)
	K_o = 5.4;   // millimolar (in standard_ionic_concentrations)
	Na_o = 140.0;   // millimolar (in standard_ionic_concentrations)
	g_K1 = 0.09;   // nanoS_per_picoF (in time_independent_potassium_current)
	tau_tr = 180.0;   // millisecond (in transfer_current_from_NSR_to_JSR)
	K_Q10 = 3.0;   // dimensionless (in transient_outward_K_current)
	g_to = 0.1652;   // nanoS_per_picoF (in transient_outward_K_current)
	phi = 0.0;

	//---------------------------------------------------------------------------
	// Computed variables
	//---------------------------------------------------------------------------

	V_rel = 0.0048 * V_cell;
	tau_u = 8.0;
	tau_f_Ca = 2.0;
	V_i = V_cell * 0.68;
	V_up = 0.0552 * V_cell;
	sigma = 1.0 / 7.0 * (exp(Na_o / 67.3) - 1.0);

	// time: time (millisecond)

	Ca_CMDN = CMDN_max * state[10] / (state[10] + Km_CMDN);
	Ca_TRPN = TRPN_max * state[10] / (state[10] + Km_TRPN);
	Ca_CSQN = CSQN_max * state[11] / (state[11] + Km_CSQN);
	i_up_leak = I_up_max * state[12] / Ca_up_max;
	i_rel = K_rel * pow(state[1], 2.0) * state[2] * state[3] * (state[11] - state[10]);
	i_Ca_L = data.GCaL * Cm * g_Ca_L * state[4] * state[6] * state[5] * (state[0] - 65.0);
	i_NaCa = Cm * I_NaCa_max * (exp(gamma * F * state[0] / (R * T)) * pow(state[14], 3.0) * Ca_o - exp((gamma - 1.0) * F * state[0] / (R * T)) * pow(Na_o, 3.0) * state[10]) / ((pow(K_mNa, 3.0) + pow(Na_o, 3.0)) * (K_mCa + Ca_o) * (1.0 + K_sat * exp((gamma - 1.0) * state[0] * F / (R * T))));
	Fn = 1.0e3 * (1.0e-15 * V_rel * i_rel - 1.0e-15 / (2.0 * F) * (0.5 * i_Ca_L - 0.2 * i_NaCa));
	u_infinity = pow(1.0 + exp(-(Fn - 3.4175e-13) / 13.67e-16), -1.0);
	// dstate[1] = (u_infinity-state[1])/tau_u;
	state[1] = Euler_inf(dt, state[1], u_infinity, tau_u);


	tau_v = 1.91 + 2.09 * pow(1.0 + exp(-(Fn - 3.4175e-13) / 13.67e-16), -1.0);
	v_infinity = 1.0 - pow(1.0 + exp(-(Fn - 6.835e-14) / 13.67e-16), -1.0);
	// dstate[2] = (v_infinity-state[2])/tau_v;
	state[2] = Euler_inf(dt, state[2], v_infinity, tau_v);

	if (fabs(state[0] - 7.9) < 1.0e-10)
		tau_w = 6.0 * 0.2 / 1.3;
	else
		tau_w = 6.0 * (1.0 - exp(-(state[0] - 7.9) / 5.0)) / ((1.0 + 0.3 * exp(-(state[0] - 7.9) / 5.0)) * 1.0 * (state[0] - 7.9));

	w_infinity = 1.0 - pow(1.0 + exp(-(state[0] - 40.0) / 17.0), -1.0);
	// dstate[3] = (w_infinity-state[3])/tau_w;
	state[3] = Euler_inf(dt, state[3], w_infinity, tau_w);
	i_up = I_up_max / (1.0 + K_up / state[10]);
	d_infinity = pow(1.0 + exp((state[0] + 10.0) / -8.0), -1.0);

	if (fabs(state[0] + 10.0) < 1.0e-10)
		tau_d = 4.579 / (1.0 + exp((state[0] + 10.0) / -6.24));
	else
		tau_d = (1.0 - exp((state[0] + 10.0) / -6.24)) / (0.035 * (state[0] + 10.0) * (1.0 + exp((state[0] + 10.0) / -6.24)));

	// dstate[4] = (d_infinity-state[4])/tau_d;
	state[4] = Euler_inf(dt, state[4], d_infinity, tau_d);
	f_Ca_infinity = pow(1.0 + state[10] / 0.00035, -1.0);
	// dstate[5] = (f_Ca_infinity-state[5])/tau_f_Ca;
	state[5] = Euler_inf(dt, state[5], f_Ca_infinity, tau_f_Ca);


	f_infinity = exp(-(state[0] + 28.0) / 6.9) / (1.0 + exp(-(state[0] + 28.0) / 6.9));
	tau_f = 9.0 * pow(0.0197 * exp(-pow(0.0337, 2.0) * pow(state[0] + 10.0, 2.0)) + 0.02, -1.0);
	// dstate[6] = (f_infinity-state[6])/tau_f;
	state[6] = Euler_inf(dt, state[6], f_infinity, tau_f);

	E_Ca = R * T / (2.0 * F) * log(Ca_o / state[10]);
	E_Na = R * T / F * log(Na_o / state[14]);
	i_B_Na = Cm * g_B_Na * (state[0] - E_Na);
	i_B_Ca = Cm * g_B_Ca * (state[0] - E_Ca);
	E_K = R * T / F * log(K_o / state[13]);
	i_B_K = Cm * g_B_K * (state[0] - E_K);



	// i_Na = Cm * g_Na * pow(state[9], 3.0) * state[7] * state[8] * (state[0] - E_Na);
	double Af = (25.2 / (1 + exp((state[0] - (-60.3+8.0)) / (-6.08))) + 16.52 / (1 + exp((state[0] - (-119.26+8.0)) / (6.21))) + 65.9) / 100.0;
	// i_Na = 3.8*Cm * g_Na * state[9] * state[9] * state[9] * (Af * state[7]*state[21] + (1.0 - Af) * state[8]*state[22]) * (state[0] - E_Na);

	double tmp = state[9]*state[9]*state[9];
	// if (tmp < 0.1)
	// {
	// 	tmp= 0.0;
	// }

	// i_Na = 1.*Cm * g_Na * tmp  * (Af * (state[7]) + (1.0 - Af) * (state[8])) * (state[0] - E_Na);
	// i_Na = 0.8*Cm * g_Na * tmp  * (Af * (tmp*state[7] + (1-tmp)*state[21]) + (1.0 - Af) * (tmp*state[8]+ (1-tmp)*state[22])) * (state[0] - E_Na);
	i_Na = 1.2 * Cm * g_Na * tmp  * (Af * (state[21]*state[7]) + (1.0 - Af) * (state[22]*state[8])) * (state[0] - E_Na);
	// i_Na = 1.5*Cm * g_Na * state[9] * state[9] * state[9] * (state[7]* state[8]) * (state[0] - E_Na);

	double V = state[0];

	double bar_m = 1.0 / (1 + exp(-(V + 46.17) / 7.0));
	double smx = (V + 45.63) / 18.81;
	double tm = (0.0855 * exp(-smx * smx)) + 0.025; //in ms

	state[9] = Euler_inf(dt, state[9], bar_m, tm);

	// double bar_h = 1. / (1 + exp((V + 79.11) / 7.02)) * 1. / (1 + exp((V + 79.11) / 7.02));
	double bar_h = 1. / (1 + exp((V + 75-5) / 6));// *1. / (1 + exp((V + 75-5) / 6));
	// double bar_h = 1. / (1 + exp((V + 75-5) / 7));
	// double bar_h = 1. / (1 + exp((V + 85) / 6.2));

	// double th1 =  30.0 / (exp((V + 66.25) / 9) * 5.0 + exp((V + 90.25) / -10.73) )  + 0.13 + 0.08;
	// double th2 = 42.838 / (exp((V + 52.295) / 11.709) * 4.29 + exp((V + 102.998) / -15.099)) + 0.43+0.3;
	// double th1 = 18.0/(exp((V+22.0)/19.5 + ((V+22.0)/32.0)*((V+22.0)/32.0)*((V+22.0)/32.0) )*60.0+ exp((V+73.0)/-18.0 ))+0.2;
	// double th2 = 62/(exp((V+6)/30 + ((V+6)/30)*((V+6)/30)*((V+6)/30) )*100.0 +exp((V+59.55)/-31.664)) +0.5;
	double Input_params [] = {
   -50.7017,
   -5.5178,
  -80.7795,
    7.3939,
   22.2304,
   33.6734,
   36.5620,
   35.3290,
    0.9774,
    0.0381,
   97.0782,
  -14.1480,
   12.9589,
   29.4515,
   80.9303,
  101.0535,
  -31.2638,
   38.3175
		};//

		double th1 = 1.0 / (exp((V + Input_params[4]) / Input_params[5] + ((V + Input_params[4]) / Input_params[5]) * ((V + Input_params[4]) / Input_params[5]) * ((V + Input_params[4]) / Input_params[5]) ) * Input_params[8] +
Input_params[9] * exp((V + Input_params[10]) / Input_params[11])) + 0.6;


double th2 = 2.82 * 62 / (exp((V + Input_params[12]) / Input_params[13] + ((V + Input_params[12]) / Input_params[13]) * ((V + Input_params[12]) / Input_params[13]) * ((V + Input_params[12]) / Input_params[13])   ) * Input_params[14]
		+ exp((V + Input_params[15]) / Input_params[16])) + 0.6 * 2.82;
// double tmp = m * m * m;

double V_shift  = V - Input_params[17];

double th1_s = 1.0 / (exp((V_shift + Input_params[4]) / Input_params[5] + ((V_shift + Input_params[4]) / Input_params[5]) * ((V_shift + Input_params[4]) / Input_params[5]) * ((V_shift + Input_params[4]) / Input_params[5]) ) * Input_params[8] +
Input_params[9] * exp((V_shift + Input_params[10]) / Input_params[11])) + 0.6;
double th2_s = 2.82 * 62 / (exp((V_shift + Input_params[12]) / Input_params[13] + ((V_shift + Input_params[12]) / Input_params[13]) * ((V_shift + Input_params[12]) / Input_params[13]) * ((V_shift + Input_params[12]) / Input_params[13])   ) * Input_params[14]
			+ exp((V_shift + Input_params[15]) / Input_params[16])) + 0.6 * 2.82;

	// double th1 = 1.0/(exp((V+26.06)/30.06 + ((V+26.06)/30.06)*((V+26.06)/30.06)*((V+26.06)/30.06) )*1.05+ 0.043*exp((V+96.221)/-18.3))+0.65;
	// double th2 =  2.82 * 62 / (exp((V + 6 + 8) / 20 + ((V + 6 + 8) / 30) * ((V + 6 + 8) / 30) * ((V + 6 + 8) / 30) ) * 100 + exp((V + 65.55 + 15) / -20.664)) + 0.6 * 2.82;

	th1 /= 2.82;
	th2 /= 2.82;
	th1_s /= 2.82;
	th2_s /= 2.82;
		// th1 = th1*4.17;

	if (state[0] < -40.0)
		alpha_h = 0.135 * exp((state[0] + 80.0) / -6.8);
	else
		alpha_h = 0.0;

	if (state[0] < -40.0)
		beta_h = 3.56 * exp(0.079 * state[0]) + 3.1e5 * exp(0.35 * state[0]);
	else
		beta_h = 1.0 / (0.13 * (1.0 + exp((state[0] + 10.66) / -11.1)));

	if (state[0] < -40.0)
		alpha_j = (-1.2714e5 * exp(0.2444 * state[0]) - 3.474e-5 * exp(-0.04391 * state[0])) * (state[0] + 37.78) / (1.0 + exp(0.311 * (state[0] + 79.23)));
	else
		alpha_j = 0.0;

	if (state[0] < -40.0)
		beta_j = 0.1212 * exp(-0.01052 * state[0]) / (1.0 + exp(-0.1378 * (state[0] + 40.14)));
	else
		beta_j = 0.3 * exp(-2.535e-7 * state[0]) / (1.0 + exp(-0.1 * (state[0] + 32.0)));

	// th2 = 1.0 / (alpha_j + beta_j);
	double ths1 =  1.0 / (alpha_h + beta_h);
	double ths2 =  1.0 / (alpha_j + beta_j);

	// state[7]+=dt*( tmp*( (1-state[7])*bar_h/th1 - state[7] * (1-bar_h)/th1  ) + (1-tmp)*  ( (1-state[7])*bar_h/ths1 - state[7] * (1-bar_h)/ths1  ) );
	// state[8]+=dt*( tmp*( (1-state[8])*bar_h/th2 - state[8] * (1-bar_h)/th2  ) + (1-tmp)*  ( (1-state[8])*bar_h/ths2 - state[8] * (1-bar_h)/ths2  ) );


	state[7] = Euler_inf(dt, state[7], bar_h, th1);
	state[8] = Euler_inf(dt, state[8], bar_h, th2);
	state[21] = Euler_inf(dt, state[21], bar_h, th1_s);//ths1*3);
	state[22] = Euler_inf(dt, state[22], bar_h, th2_s);


	/*	if (state[0] < -40.0)
			alpha_h = 0.135 * exp((state[0] + 80.0) / -6.8);
		else
			alpha_h = 0.0;

		if (state[0] < -40.0)
			beta_h = 3.56 * exp(0.079 * state[0]) + 3.1e5 * exp(0.35 * state[0]);
		else
			beta_h = 1.0 / (0.13 * (1.0 + exp((state[0] + 10.66) / -11.1)));

		h_inf = alpha_h / (alpha_h + beta_h);
		tau_h = 1.0 / (alpha_h + beta_h);
		// dstate[7] = (h_inf-state[7])/tau_h;
		state[7] = Euler_inf(dt, state[7], h_inf, tau_h);

		if (state[0] < -40.0)
			alpha_j = (-1.2714e5 * exp(0.2444 * state[0]) - 3.474e-5 * exp(-0.04391 * state[0])) * (state[0] + 37.78) / (1.0 + exp(0.311 * (state[0] + 79.23)));
		else
			alpha_j = 0.0;

		if (state[0] < -40.0)
			beta_j = 0.1212 * exp(-0.01052 * state[0]) / (1.0 + exp(-0.1378 * (state[0] + 40.14)));
		else
			beta_j = 0.3 * exp(-2.535e-7 * state[0]) / (1.0 + exp(-0.1 * (state[0] + 32.0)));

		j_inf = alpha_j / (alpha_j + beta_j);
		tau_j = 1.0 / (alpha_j + beta_j);
		// dstate[8] = (j_inf-state[8])/tau_j;
		state[8] = Euler_inf(dt, state[8], j_inf, tau_j);

		if (state[0] == -47.13)
			alpha_m = 3.2;
		else
			alpha_m = 0.32 * (state[0] + 47.13) / (1.0 - exp(-0.1 * (state[0] + 47.13)));

		beta_m = 0.08 * exp(-state[0] / 11.0);
		m_inf = alpha_m / (alpha_m + beta_m);
		tau_m = 1.0 / (alpha_m + beta_m);
		// dstate[9] = (m_inf-state[9])/tau_m;
		state[9] = Euler_inf(dt, state[9], m_inf, tau_m);
	*/

	f_NaK = pow(1.0 + 0.1245 * exp(-0.1 * F * state[0] / (R * T)) + 0.0365 * sigma * exp(-F * state[0] / (R * T)), -1.0);
	i_NaK = Cm * i_NaK_max * f_NaK * 1.0 / (1.0 + pow(Km_Na_i / state[14], 1.5)) * K_o / (K_o + Km_K_o);
	// dstate[14] = (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F);
	state[14] += dt * ((-3.0 * i_NaK - (3.0 * i_NaCa + i_B_Na + i_Na)) / (V_i * F));

	i_K1 = data.GK1 * Cm * g_K1 * (state[0] - E_K) / (1.0 + exp(0.07 * (state[0] + 80.0)));
	i_to = data.Gto * Cm * g_to * pow(state[17], 3.0) * state[18] * (state[0] - E_K);
	g_Kur = 0.005 + 0.05 / (1.0 + exp((state[0] - 15.0) / -13.0));
	i_Kur = data.Simple_GKur * Cm * g_Kur * pow(state[19], 3.0) * state[20] * (state[0] - E_K);
	i_Kr = data.GKr * Cm * g_Kr * state[15] * (state[0] - E_K) / (1.0 + exp((state[0] + 15.0) / 22.4));
	i_Ks = data.GKs * Cm * g_Ks * pow(state[16], 2.0) * (state[0] - E_K) + phi * Cm * g_Ks * pow(state[16], 2.0) * (state[0] - E_K);
	// dstate[13] = (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F);
	state[13] += dt * ((2.0 * i_NaK - (i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_K)) / (V_i * F));


	i_CaP = Cm * i_CaP_max * state[10] / (0.0005 + state[10]);
	B1 = (2.0 * i_NaCa - (i_CaP + i_Ca_L + i_B_Ca)) / (2.0 * V_i * F) + (V_up * (i_up_leak - i_up) + i_rel * V_rel) / V_i;
	B2 = 1.0 + TRPN_max * Km_TRPN / pow(state[10] + Km_TRPN, 2.0) + CMDN_max * Km_CMDN / pow(state[10] + Km_CMDN, 2.0);
	// dstate[10] = B1/B2;
	state[10] += dt * B1 / B2;

	i_tr = (state[12] - state[11]) / tau_tr;
	// dstate[12] = i_up-(i_up_leak+i_tr*V_rel/V_up);
	state[12] += dt * (i_up - (i_up_leak + i_tr * V_rel / V_up));
	// dstate[11] = (i_tr-i_rel)*pow(1.0+CSQN_max*Km_CSQN/pow(state[11]+Km_CSQN, 2.0), -1.0);
	state[11] += dt * ((i_tr - i_rel) * pow(1.0 + CSQN_max * Km_CSQN / pow(state[11] + Km_CSQN, 2.0), -1.0));

	if (fabs(state[0] + 14.1) < 1.0e-10)
		alpha_xr = 0.0015;
	else
		alpha_xr = 0.0003 * (state[0] + 14.1) / (1.0 - exp((state[0] + 14.1) / -5.0));

	if (fabs(state[0] - 3.3328) < 1.0e-10)
		beta_xr = 3.7836118e-4;
	else
		beta_xr = 0.000073898 * (state[0] - 3.3328) / (exp((state[0] - 3.3328) / 5.1237) - 1.0);

	tau_xr = pow(alpha_xr + beta_xr, -1.0);
	xr_infinity = pow(1.0 + exp((state[0] + 14.1) / -6.5), -1.0);
	// dstate[15] = (xr_infinity-state[15])/tau_xr;
	state[15] = Euler_inf(dt, state[15], xr_infinity, tau_xr);

	if (fabs(state[0] - 19.9) < 1.0e-10)
		alpha_xs = 0.00068;
	else
		alpha_xs = 0.00004 * (state[0] - 19.9) / (1.0 - exp((state[0] - 19.9) / -17.0));

	if (fabs(state[0] - 19.9) < 1.0e-10)
		beta_xs = 0.000315;
	else
		beta_xs = 0.000035 * (state[0] - 19.9) / (exp((state[0] - 19.9) / 9.0) - 1.0);

	tau_xs = 0.5 * pow(alpha_xs + beta_xs, -1.0);
	xs_infinity = pow(1.0 + exp((state[0] - 19.9) / -12.7), -0.5);
	// dstate[16] = (xs_infinity-state[16])/tau_xs;
	state[16] = Euler_inf(dt, state[16], xs_infinity, tau_xs);



	alpha_oa = 0.65 * pow(exp((state[0] - (-10.0)) / -8.5) + exp((state[0] - (-10.0) - 40.0) / -59.0), -1.0);
	beta_oa = 0.65 * pow(2.5 + exp((state[0] - (-10.0) + 72.0) / 17.0), -1.0);
	tau_oa = pow(alpha_oa + beta_oa, -1.0) / K_Q10;
	oa_infinity = pow(1.0 + exp((state[0] - (-10.0) + 10.47) / -17.54), -1.0);
	// dstate[17] = (oa_infinity-state[17])/tau_oa;
	state[17] = Euler_inf(dt, state[17], oa_infinity, tau_oa);




	alpha_oi = pow(18.53 + 1.0 * exp((state[0] - (-10.0) + 103.7) / 10.95), -1.0);
	beta_oi = pow(35.56 + 1.0 * exp((state[0] - (-10.0) - 8.74) / -7.44), -1.0);
	tau_oi = pow(alpha_oi + beta_oi, -1.0) / K_Q10;
	oi_infinity = pow(1.0 + exp((state[0] - (-10.0) + 33.1) / 5.3), -1.0);
	// dstate[18] = (oi_infinity-state[18])/tau_oi;
	state[18] = Euler_inf(dt, state[18], oi_infinity, tau_oi);

	alpha_ua = 0.65 * pow(exp((state[0] - (-10.0)) / -8.5) + exp((state[0] - (-10.0) - 40.0) / -59.0), -1.0);
	beta_ua = 0.65 * pow(2.5 + exp((state[0] - (-10.0) + 72.0) / 17.0), -1.0);
	tau_ua = pow(alpha_ua + beta_ua, -1.0) / K_Q10;
	ua_infinity = pow(1.0 + exp((state[0] - (-10.0) + 20.3 - data.Simple_Ikur_ac_shift) / (-9.6 * data.Simple_Ikur_ac_grad)), -1.0);
	// dstate[19] = (ua_infinity-state[19])/tau_ua;
	state[19] = Euler_inf(dt, state[19], ua_infinity, tau_ua);

	alpha_ui = pow(21.0 + 1.0 * exp((state[0] - (-10.0) - 195.0) / -28.0), -1.0);
	beta_ui = 1.0 / exp((state[0] - (-10.0) - 168.0) / -16.0);
	tau_ui = pow(alpha_ui + beta_ui, -1.0) / K_Q10;
	ui_infinity = pow(1.0 + exp((state[0] - (-10.0) - 109.45 - data.Simple_Ikur_inac_shift) / (27.48 * data.Simple_Ikur_inac_grad)), -1.0);
	// dstate[20] = (ui_infinity-state[20])/tau_ui;
	state[20] = Euler_inf(dt, state[20], ui_infinity, tau_ui);

	// state[22] = i_Na ;


	/* New IKur for CNZ model */
	/*double const    CNZ_gkur = 0.006398;
	double V = state[0];
	double CNZ_a = state[19];
	double CNZ_i = state[20];
	//IKur
	i_Kur = data.GKur * Cm * CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597)))) * CNZ_a * CNZ_i * (V - E_K);
	double K_Q10_CNZ = 3.5308257834747638;
	//CNZ_a
	double inf = ((data.IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + data.IKur_ac1_shift) / (-5.75418 * data.IKur_ac1_grad))) ) * ((data.IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + data.IKur_ac2_shift) / (-11.51037561 * data.IKur_ac2_grad)))) + data.IKur_ac_add;
	double tau = (45.6666746826 / (1 + exp((V + 11.2306497073) / 11.5254705962)) + 4.26753514993) * (0.262186042981 / (1 + exp((V + 35.8658312707) / (-3.87510627762))) + 0.291755017928); //
	tau = tau / K_Q10_CNZ;

	// CNZ_a = inf + (CNZ_a - inf) * exp(-(dt) / tau);
	// dY[19] = (inf - Y[19])/tau;
	state[19] = Euler_inf(dt, state[19], inf, tau);

	// CNZ_i
	inf = (data.IKur_inac_mult * 0.52424) / (1.0 + exp((V + 15.1142 + data.IKur_inac_shift2) / (7.567021 * data.IKur_inac_grad2))) + 0.4580778 + data.IKur_inac_add;
	tau = 2328 / (1 + exp(((V) - 9.435) / (3.5827))) + 1739.139;
	tau = tau / K_Q10_CNZ;
	// CNZ_i = inf + (CNZ_i - inf) * exp(-(dt) / tau);
	// dY[20] = (inf - Y[20])/tau;
	state[20] = Euler_inf(dt, state[20], inf, tau);*/
	// for debug purpose...
	/*for (int i = 0; i < 21; ++i)
	{
		if (state[i] != state[i])
		{
			std::cerr << i << " nan" << std::endl;
			std::exit(0);
		}
	}
	std::cerr << state[0] << ' ' << state[8] << std::endl;*/
	return (i_Na + i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP + i_NaCa + i_Ca_L ) / Cm + Istim;
}
