/*
 * new_atrium_LAA_RAA.c
 * new atrium code with LAA and RAA included, geometry has no AVN, SAN modelled as SAN or CT
 *
 */
#include <memory>
#include <algorithm>
#include <vector>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// #include "util.h"
// #include "conduction.h"
// #include "tissue.h"
// #include "cell.h"
#include "stimulus.h"
#include "input_output.h"
#include <string.h>
#include "simulation_config.h"

#include <exception>

// #include "Tent_Vent.h"
// #include "CNZ.h"
// #include "Updated_CNZ.h"

#include "NCZ_ODE.h"

// #include "rudy.h"
#include "EnumSimulationCtrl.hpp"

// #include "TNNP_MarkovIKr_SP.hpp"
// #include "TNNP_MarkovIKr_function.hpp"
// #include "TNNP_ISO.hpp"
// #include "TNNP_MarkovIKr.h"
#include "one_d_conduction.h"
#include "APInfo.hpp"

// if we're compiling with openMP, include the file, else define some stub
#include <omp.h>
// #define _OPENMP 1

// functions
#ifdef _OPENMP
#include <omp.h>
#else
// #define omp_get_num_threads() (48)
// #define omp_get_thread_num() (0)
#endif
// #define NX (235)
// #define NY (269)
// #define NZ (298)

#define OUT_TO_TXT (1)
// #define BREAK_IF_CVM // break if CV measured values reach steady-state values.
// #define OUT_INICOND
/* define the num of neighboorhood ( related to the model dimension */
#define CV_COM (1)
#define NUM_NBH (3)
#define DIFF_COEF (0.15*1.525*0.63*1.45*(0.6))  // *0.6 40% reduction // wild type, CV = 0.752688


// #define DIFF_COEF (0)
#define FIRST_ACT (100)
#define SECOND_ACT (NX-20)

/* Downsampled IS ventricle dimension */
#define NX (200) // 200 cells in a single strand
#define NY (325)
#define NZ (358)

// #define PROFILE_TIME
/*#define NX (100)
#define NY (100)
#define NZ (100)*/


// #define D2 (0.21)
// #define DD (8*0.21)
/* to be determined later */
/* The diffusion parameter */
#define D2 (0.18/10.0)
#define DD (0.18)

// typedef TTCell TNNP_MarkovIKr;

#define STATE_LEN (40)//cnz_con_length())  /* the state length of Tent model is 21,crn_con_length() = 39*/


int main (int argc, char **argv) {
	omp_set_num_threads(8);
	// omp_set_dynamic(1);
	// counters
	int i, count = 0, c;

	// tissue variables
	unsigned char           *** atrium;
	unsigned char            *tissue;

	// propagation variables
	const double            dx = 0.30;
	int                      *nbd;
	int                      **neighbourhood;
	double                   *lap;
	double                   **laplacian;
	double                   *dv;
	double                   *v_new;
	double                   *v_old;
	double 					 *dv_diff_temp;
	double 					 *RK_K1;
	double 					 *RK_K2;
	double 					 *RK_K3;
	double 					 *v_temp;
	double 					 *v_temp_2;


	const int S1_number = 10;

	// time variables
	const double dt      = 0.1;  // large Dt here, addpative dt for single cell computations.
	double  t_max        = 2000.0;
	double  t;

	Simulation_Config Config_para(argc, argv);
	Config_para.dt = dt;
	Config_para.Report_Config();
// #ifdef OUT_INICOND
// 	const double t_start = 0.0;//(Config_para.S1_number-1)*Config_para.BCL - 10;
// #else
// 	const double t_start = - 10 + dt;//(Config_para.S1_number - 1) * Config_para.BCL - 10 + dt;
// #endif
	// cell stuff
	double *state;
	// lookup_table_t           *table;
	double                  results[1];
	int conduction = 1;
	// std::vector<TNNP_MarkovIKr> cell_vec;

	// stim stuff
	const double stim   = 20.0;
	char         *san   = NULL;
	char         *sanin = NULL;
	float        *stims = NULL;
	gzFile       gz     = NULL;
	int          iStim  = 0;

	// output variables
	int       outcount  = 0;
	int       filecount = 0;
	int       outindex  = 0;
	float     *v_out;

	float *Fcell;
	FILE *in;
	int *FB_number;
	int FB_type;

	float ISO, ACH;
	int AF;

	int BCL, S2;
	double stimint;
	int intBCL, intstim, stimcount, timeint, durationint;
	double stimduration;
	int stimflag;
	double Dscale    = 1.0;

	double Ggap;

	int IKur_switch  = 0;

	// #ifdef _OPENMP
	double start     = omp_get_wtime();
	double end       = 0;
	double end_setup = 0;
	// #endif
	float *Stim_amp, * Stim_time, *Stim_amp_S2;

	std::vector<int> AP_Number(NX);
	// Dealing with Config_paraments etc
	// Default values

	// Config_para.Total_time = (Config_para.S1_number) * Config_para.BCL + Config_para.S2 * 2 + 100.0;
	std::vector<double> PrePotential(NX);
	std::vector<double> TemPotential(NX);
	std::vector<double> ActTime(NX);
	ISO     = ACH = AF = 0;
	BCL     = 1000;
	S2      = 0;
	FB_type = 0;
	Ggap    = 3.0;
	double timeelasped;
	SingleCellPara *cell_para;
	std::ofstream OneD_output;
	std::ofstream cv_out_file("cv_out_file.dat", std::ios::out);
	std::vector<TypeCell> cell_type_vec;
	std::vector<double> cv; // for recording cv
	std::vector<std::shared_ptr<APInfor> > APInfor_vec;

	// omp_set_num_threads(4);


	for (auto& x : Config_para.Popul_scalingF)
		std::cout << x << ' ';



	// assign scaling factors to conductances.


	double t_first_act, t_second_act;
	#pragma omp parallel  default (none) \
	shared(atrium, stderr, Config_para, \
	       v_new, v_old,dv_diff_temp, v_temp, v_temp_2, RK_K3,RK_K2,RK_K1, std::cerr, \
	       neighbourhood, tissue, laplacian, c,conduction, end_setup, sanin, Fcell, gz, in,  cell_type_vec, \
	       argv, ISO, ACH, AF, argc, intstim, stimint, intBCL, \
	       durationint, iStim, S2, BCL, FB_number, FB_type, OneD_output,\
	       Ggap, IKur_switch, Dscale, t_max, start,AP_Number,PrePotential, TemPotential, end,timeelasped,\
	       Stim_time, Stim_amp,Stim_amp_S2, stims, san,state, t_first_act, ActTime, t_second_act,cv,cell_para,APInfor_vec,nbd, lap) \
	private( dv,  t, i, count, v_out,  \
	         results, stimcount, stimflag, timeint) \
	firstprivate(outcount, filecount, outindex)
	{

		#pragma omp single
		{
			/*3D*/
			/*atrium = read_and_embed_geometry(Config_para.Geometry_file, NZ, NY, NX);
			neighbourhood = generate_neighbours_map(atrium, NZ, NY, NX, &c);
			laplacian = generate_laplacian(atrium, NZ, NY, NX, c,
			Config_para.Fibre_theta, Config_para.Fibre_phi, dx, DD, D2);
			tissue = generate_tissue_map(atrium, NZ, NY, NX, c);*/
			/*laplacian = generate_laplacian_test_heterogeneity(atrium, NZ, NY, NX, c,
			Fibre_theta, Fibre_phi, dx, DD, D2);*/

			/*laplacian = generate_laplacian_using_fibre_components(atrium, NZ, NY, NX, c,
			"Last_Atria_Geo_Full_Fibre_X.geo.gz", "Last_Atria_Geo_Full_Fibre_Y.geo.gz", "Last_Atria_Geo_Full_Fibre_Z.geo.gz",
			dx, DD, D2);*/
			/* for One Dimensional */
			tissue = new unsigned char [NX];
			for (int id = 0; id < NX; id++) {
				tissue[id] = 14;
			}
			c = NX;
			std::cerr << ">>Succesfully read in Geometry and Fibre Files<<\n";
			std::cerr << ">>Total Number of Cells: " << c << " <<\n";
			laplacian = generate_one_D_laplacian_test(NX, tissue, DIFF_COEF, dx);
			neighbourhood = generate_one_D_neighboors(NX);
			// OneD_output = fopen("OneD_output.dat", "wt");
			OneD_output.open("OneD_output.dat");
			//output_matrix(laplacian, c, 19);
		}




		// create local copies of the neighbourhood and laplacian.
		const int count = c;


		// MALLOC(v_out, count * sizeof(float));
		v_out = new float [count];
		// dv    = new double [count];
		// MALLOC(dv, count * sizeof(double));


		// allocate voltage arrays
		#pragma omp single
		{



			nbd = new int [c * NUM_NBH];
			lap = new double [c * NUM_NBH];


			for (int n = 0; n < count; n++) {
				for (i = 0; i < NUM_NBH; i++) {
					nbd[(n * NUM_NBH) + i] = neighbourhood[n][i];
					lap[(n * NUM_NBH) + i] = laplacian[n][i];
				}
			}


			state = new double [ STATE_LEN * count];
			v_new     = new double [count];
			v_old     = new double [count];
			dv_diff_temp = new double [count];
			v_temp = new double [count];
			v_temp_2 = new double [count];
			RK_K1 = new double [count];
			RK_K2 = new double [count];
			RK_K3 = new double [count];
			san       = new  char [count];
			stims     = new float [count];
			FB_number = new int [count];
			Fcell     = new float [count];
			Stim_time = new float [count];
			Stim_amp  = new float [count];
			Stim_amp_S2 = new float[count];
			cell_para = new SingleCellPara [count];
			// cell_vec  = new TNNP_MarkovIKr* [count];
			/*for (int i = 0; i < count; ++i)
			{
			    TypeCell celltype = EPI;
			    double epi_mid_ratio = 1.5;
			    TNNP_MarkovIKr cell( celltype, epi_mid_ratio, dt );
			    cell_vec.push_back(cell);
			}*/

			/*read_float_from_bin(Config_para.Stim_time_file, Stim_time, count);
			read_float_from_bin(Config_para.Stim_amp_file, Stim_amp, count);*/

			for (int i = 0; i < count; ++i) {
				Stim_time[i] = 0.0;
				Stim_amp[i]  = 20.0;
				Stim_amp_S2[i] = 0.0;
				if (i > 15)
					Stim_amp[i] = 0.0;

				// if (i > 102 and i < 113)
				if (i < 16)
					Stim_amp_S2[i] = 30.0;
			}

			if (!Config_para.SAN_type.c_str())
				perror("SAN Setup wrong");
			/*fprintf(stderr, "SAN_type = %s\n", Config_para.SAN_type);
			fprintf(stderr, "Stim_time_File = %s\n", Config_para.Stim_time_file);*/
			std::cerr << "SAN_type = " << Config_para.SAN_type << std::endl;
			std::cerr << "Stim_time_File = " << Config_para.Stim_time_file << std::endl;
			// printf(">>Assigned Fcell according to parameters<<\n");
			std::cerr << ">>Assigned Fcell according to parameters<<" << std::endl;


			// FB stuff
			if (FB_type != 0) {
				in = fopen(Config_para.FB_map.c_str(), "r");
				for (int n = 0; n < count; n++) {
					// fscanf(in, "%d\n", &FB_number[n]);
				}
				std::cerr << ">>Assigned FB map<<<<" << std::endl;
			} else {
				for (int n = 0; n < count; n++) FB_number[n] = 0;
			}


			for (int i = 0; i < count; ++i)
			{
				AP_Number[i] = 0; // initial beats numbers.

				if (i == SECOND_ACT)
				{
					APInfor_vec.push_back(std::shared_ptr<APInfor>(new APInfor("APD2.txt", false, true) ));
				} else if (i == FIRST_ACT)
				{
					APInfor_vec.push_back(std::shared_ptr<APInfor>(new APInfor("APD1.txt", false, true) ));
				}

				else {
					// APInfor temp;
					APInfor_vec.push_back(std::shared_ptr<APInfor>(new APInfor() ));
				}
			}
		} // end omp single


		// assign voltages
		// #pragma omp for schedule(static)
		#pragma omp single
		for (int n = 0; n < count; n++) {
			// Regionalparameters(&cell_para[n], 1e-15, 1e-15, tissue[n]);

			IKurMuationType mutation = WT;
			if (Config_para.mutation_char == "WT" or  Config_para.mutation_char == "None") {
				mutation = WT;
			} else if (Config_para.mutation_char == "D322H" ) {
				mutation = D322H;
			} else if (Config_para.mutation_char == "E48G" ) {
				mutation = E48G;
			} else if (Config_para.mutation_char == "A305T" ) {
				mutation = A305T;
			}

			else if (Config_para.mutation_char == "Y155C" ) {
				mutation = Y155C;
			} else if (Config_para.mutation_char == "D469E" ) {
				mutation = D469E;
			} else if (Config_para.mutation_char == "P488S" ) {
				mutation = P488S;
			}

			else {
				std::cerr << "mutation wrong!!!!\n";
				std::exit(0);
			}

			// AtriaMutationParameters(&cell_para[n], mutation);
			TypeCell type;
			type = GetCellTypeFromLabel(14);//tissue[n]);

			cell_para[n].SetCellTypePara(type);
			cell_para[n].SetIKurMutationPara(mutation);
			// cell_para[n].drug_INa_concen = 0.0;

			// cell_para[n].SetINaBlockPara(60 * 1e-6, 1, 200, 0.760035 , 0.007203);

			/* assign drug effects here if needed */
			// cell_para[n].SetINaBlockPara(60 * 1e-6, 100.0, 100.0, 1, 0.01);

			// cell_para[n].SetAcacetinEffectPara(3.2); // 3.2 um of acacetin
			// cell_para[n].SetAcacetin_IKurOnly(3.2);

			if (Config_para.AF_model > 0) {
				cell_para[n].SetAFTypePara(AF4);
			}

			if (Config_para.Popl_SF_Number == 15) {
				cell_para[n].GNa   *= Config_para.Popul_scalingF[0];
				cell_para[n].GNaL  *= Config_para.Popul_scalingF[1];
				cell_para[n].GCaL  *= Config_para.Popul_scalingF[2];
				cell_para[n].Gto   *= Config_para.Popul_scalingF[3];
				cell_para[n].GKur  *= Config_para.Popul_scalingF[4];
				cell_para[n].GKr   *= Config_para.Popul_scalingF[5];
				cell_para[n].GKs   *= Config_para.Popul_scalingF[6];
				cell_para[n].GK1   *= Config_para.Popul_scalingF[7];
				cell_para[n].GNaCa *= Config_para.Popul_scalingF[8];
				cell_para[n].GNaK  *= Config_para.Popul_scalingF[9];
				cell_para[n].GSK   *= Config_para.Popul_scalingF[10];
				cell_para[n].GK2P  *= Config_para.Popul_scalingF[11];
				cell_para[n].GbNa  *= Config_para.Popul_scalingF[12];
				cell_para[n].GbCa  *= Config_para.Popul_scalingF[13];
				cell_para[n].GCap  *= Config_para.Popul_scalingF[14];
			}

			if (Config_para.Remodelling_F_Number == 3 and (Config_para.AF_model == 4) ) {
				cell_para[n].GKur  *= Config_para.Remodelling_F[0];
				cell_para[n].GSK   *=  Config_para.Remodelling_F[1];
				cell_para[n].GK2P  *= Config_para.Remodelling_F[2];
			}

			if (Config_para.Modulation_F_Number == 3) {
				cell_para[n].GKur  *= 1.0 - Config_para.Modulation_F[0];
				cell_para[n].GSK   *=  1.0 - Config_para.Modulation_F[1];
				cell_para[n].GK2P  *= 1.0 - Config_para.Modulation_F[2];
			}



			// cell_para[n].SetAFTypePara(AF3);
			// set cell types heterogeneity
			// unsigned char label = tissue[n];
			// TypeCell celltype = LVENDO;

			// // set initial conditions
			// double epi_mid_ratio = 1.5;
			// // crn_con_initial_conditions_Vent(&(state[n * STATE_LEN]), tissue[n]);
			// // cell_vec[n] = new TNNP_MarkovIKr_SP( celltype, epi_mid_ratio, dt );
			// // cell_vec.at(n).SetCellType(celltype);
			// // crn_con_initial_conditions_Vent(&(state[n * STATE_LEN]), tissue[n]);
			// // set voltage of cell equal to its initial condition voltage
			// // v_new[n] = v_old[n] = cell_vec.at(n).GetPotential();
			// if (n < 20)
			// 	celltype = LVENDO;
			// else if (n > 20 && n < 70)
			// 	celltype = LVMCELL;
			// else
			// 	celltype = LVEPI;


			// printf("%s\n", );

			// std::cerr << "ya"
			if (Config_para.ICs == "Default") {
				/* code */
				// INa_IKur_Drug_Updated_CNZ_ODE_Initialise(&(state[n * STATE_LEN]), tissue[n]);
				NCZ_ODE_Update_Initialise(&(state[n * STATE_LEN]), tissue[n]);


			} else if (Config_para.ICs == "Specific") {

				char filename [100];
				sprintf(filename, "IniCond/Cell_%03d_BCL_%04d.bin", n, (int) Config_para.BCL);
				// sprintf(filename, "Restart_ICs/AF_CTL/Drug_BCL_%d_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.bin", int(Config_para.BCL), 14, "WT", Config_para.AF_model);
				// sprintf(filename, "Restart_ICs/AF_INa_Y/Drug_BCL_%d_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.bin", int(Config_para.BCL), 14, "WT", Config_para.AF_model);
				// sprintf(filename, "Restart_ICs/IKur_3.2_INa_NO/Drug_BCL_%d_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.bin", int(Config_para.BCL), 14, "WT", Config_para.AF_model);
				// sprintf(filename, "Restart_ICs/IKur_3.2_INa_Y/Drug_BCL_%d_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.bin", int(Config_para.BCL), 14, "WT", Config_para.AF_model);
				// sprintf(filename, "Restart_ICs/Acacetin_3.2_INa_Y/Drug_BCL_%d_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.bin", int(Config_para.BCL), 14, "WT", Config_para.AF_model);
				// sprintf(filename, "Restart_ICs/Acacetin_3.2_INa_Y/Drug_BCL_%d_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.bin", 170, 14, "WT", Config_para.AF_model);
				// std::cerr << filename;

				// sprintf(filename, "Restart_ICs/ICs.bin");  // get data from the model.
				read_double_from_bin(filename, &state[n * STATE_LEN], STATE_LEN);
			} else {
				read_double_from_bin(Config_para.ICs.c_str(), &state[n * STATE_LEN], STATE_LEN);
			}

			// cell_type_vec.push_back(celltype);
			// TNNP_MarkovIKr_Initialise(&(state[n * STATE_LEN]), celltype);
			// InitialseTentVentByJohn(&(state[n * STATE_LEN]), celltype);
			v_new[n] = v_old[n] = state[n * STATE_LEN];
#ifdef CV_COM
			PrePotential[n] = v_new[n];
			TemPotential[n] = v_new[n];
			ActTime[n] = -1000.0;
#endif
			/*regional cell models parameter assignment */
			/* to be continued.... */
			/* this part was updated in the latest wholeheart main function */
			/* if you want to complete this part, please refer to the new file system Paramters.c*/
		}


		// #ifdef _OPENMP
		#pragma omp master
		{
			end_setup = omp_get_wtime();
			// fprintf(stderr, "Setup time:   %6g\n", end_setup - start);
			std::cerr << "Setup time:  " << end_setup - start << std::endl;
		}
		// #endif
		#pragma omp barrier
		#pragma omp single
		{
			// fprintf(stderr, ">>Assigned regional and modulation parameters<<\n>>Starting time loop now<<\n");
			std::cerr << ">>Assigned regional and modulation parameters<<\n>>Starting time loop now<<" << std::endl;
		}
		#pragma omp barrier

		double output_time = 0.0;
		timeelasped = omp_get_wtime();

		// loop over time.
		for (t = Config_para.t_start; t < Config_para.Total_time; t += dt) {

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; ++n) {

#ifdef OUT_INICOND

				if (t >= Config_para.Total_time - 2 * Config_para.BCL - 10 - dt / 2.0 and t <= Config_para.Total_time - 2 * Config_para.BCL - 10 + dt / 2.0) // output initial values
				{
					char filename [100];
					sprintf(filename, "IniCond/Cell_%03d_BCL_%04d.bin", n, (int) Config_para.BCL);
					output_double_array_bin(filename, &state[n * STATE_LEN], STATE_LEN);
				}
#endif
			}





			if (!outcount) {
				// if the thread index is the current thread, store a copy of the
				// state.
				if (outindex == omp_get_thread_num()) {
#ifdef PROFILE_TIME

					std::cerr << "in time loop now, t = " << t << " <<" << std::endl;
					std::cerr << "time spent = " << omp_get_wtime() - timeelasped << " <<" << std::endl;
					timeelasped = omp_get_wtime();
#endif
					for (int n = 0; n < count; n++)
						v_out[n] = (float) v_new[n];
				}
				// increment various counters.
				outindex++;
				outcount = 1;//(int) (0.1/ dt);

#ifdef OUT_TO_TXT
				#pragma omp single
				{
					if (t > Config_para.Total_time - 10 * Config_para.BCL - 10)
					{
						OutPutArrayToTxt(OneD_output, v_new, count);
					}
				}
#endif
				#pragma omp barrier
				// if all the threads have output.
				if (outindex == omp_get_num_threads()) {
#ifdef OUT_TO_BIN
					output_voltage_array(count, v_out, "v%04d.bin", (filecount * omp_get_num_threads()) + omp_get_thread_num());
					filecount++;
#endif
					outindex = 0;
				}
			}
			outcount--;

			// stops threads from potentially updating something before it's been
			// output.

			/*#pragma omp for schedule(static)
			for (n = 0; n < count; n++) {
			    dv[n] = 0.0;
			    stims[n] = 0.0;
			    for (i = 0; i < 19; i++)
			        dv[n] += (lap[(n * 19) + i] * v_old[nbd[(n * 19) + i]]);
			}
			*/

#ifdef PROFILE_TIME
			#pragma omp single
			{
				std::cerr << "Setup time:  " <<  omp_get_wtime() - end_setup << std::endl;
				timeelasped = omp_get_wtime();
			}
#endif

			/*operator spliting part 1*/
			/*#pragma omp for schedule(static)
			for (int n = 0; n < count; n++) {
				double sum = 0.0;
				// #pragma omp simd reduction // not available on icpc 12.0
				for (int ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_old[nbd[(n * NUM_NBH) + ii]]);
				}
				v_new[n] = v_old[n] + sum * dt / 2.0;
			}*/


			/*operator spliting part 1*/

			/* using Heun's Method -> https://en.wikipedia.org/wiki/Heun%27s_method */
			/*#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_new[n] = v_old[n];  // keep a copy of v_old[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_old[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				dv_diff_temp[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_old[n] + sum * dt / 2.0;  // devide by two because of operator splitting.
			}

			#pragma omp for schedule(static)   // part 2, y_i+1 = y_i + h/2*(f(t, y(t)) + f(t+1, y_bar_i+1))
			for (int n = 0; n < count; n++) {
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;
				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]);  // this is f(t+1, y_bar_i+1)
				}
				//y_bar_i+1
				v_new[n] = v_new[n] + (dt / 2.0) / 2.0 * (dv_diff_temp[n] + sum);  // devide by two because of operator splitting.
			}
			*/

			/* using two step RK https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods*/
			/*#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_new[n] = v_old[n];  // keep a copy of v_old[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_old[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				dv_diff_temp[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_old[n] + sum * (dt / 2.0) / 2.0;  // devide by two because of operator splitting.
			}

			#pragma omp for schedule(static)   // part 2, y_i+1 = y_i + h/2*(f(t, y(t)) + f(t+1, y_bar_i+1))
			for (int n = 0; n < count; n++) {
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;
				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]);  // this is f(t+1, y_bar_i+1)
				}
				// dv[n] = sum;
				//y_bar_i+1
				v_new[n] = v_new[n] + (dt / 2.0) *( sum);  // devide by two because of operator splitting.
			}*/


			/*RK4 method*/
			/* using RK4 here*/
			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_new[n] = v_old[n];  // keep a copy of v_old[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_old[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K1[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_old[n] + sum * dt / 2.0 / 2.0; // devide by two because of operator splitting. // y_n_K1/2
			}
			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_new[n] = v_old[n];  // keep a copy of v_old[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K2[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp_2[n] = v_old[n] + sum * dt / 2.0 / 2.0; // devide by two because of operator splitting. // y_n_K2/2
			}

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_new[n] = v_old[n];  // keep a copy of v_old[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp_2[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K3[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_old[n] + sum * dt / 2.0; // devide by two because of operator splitting. // y_n_K3/2
			}

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_new[n] = v_old[n];  // keep a copy of v_old[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				v_new[n] = v_new[n] + (sum + RK_K1[n] + 2 * RK_K2[n] + 2 * RK_K3[n]) * dt / 2.0 / 6.0; // devide by two because of operator splitting. // y_n_K3/2
				// APInfor_vec[n]->MeasureAPD90Using_dVdtMax(t, 0, Config_para.BCL, dt, v_new[n]);  // measure APD at the end of each time step.
			}






			#pragma omp for schedule(static)
			for (int n = 0; n < count; n++) {
				// stims[n] = iStim * stim * san[n]; // stim = iStim (stim time yes/no) * stim size * san[n] (stim location yes/no)
				/* stimuli */
				// standard S1 stimuli
				// stims[n] = S1(1000 * Stim_time[n], Stim_amp[n], BCL, t, 2.0);
				// cell_vec.at(n).svolt = v_new[n];
				/*cell_vec.at(n).SetPotential(v_new[n]);
				cell_vec.at(n).RunSingleTimeStep(dt, stims[n]);
				cell_vec.at(n).ReturnPotentialByReference(v_new[n]);*/
				// double stims = S1S2_num(0.0, Stim_amp[n], Config_para.BCL, Config_para.S1_number, Config_para.S2, 1, t, 2.0);
				// double stims = S1S2_num(0, Stim_amp[n], Config_para.BCL, Config_para.S1_number, Config_para.S2, 1, t, 2.0);

				double stims = S1S2_num_stim(0, Config_para.BCL, Stim_amp[n], Config_para.S1_number,
				                             Config_para.S2, Stim_amp_S2[n], 1, t, 2.0);


				int K_factor = 1 + (int) (fabs(APInfor_vec[n]->dVdt_pre));
				K_factor = K_factor > 5 ? 5 : K_factor;

				double state_temp[STATE_LEN]; // make a copy of state variables
				for (int jj = 1; jj < STATE_LEN; ++jj) {
					state_temp[jj] = state[(n * STATE_LEN) + jj];   // make a copy of state variables
				}

				if (K_factor > 1) {
					double small_dt = dt / K_factor;
					for (int ii = 0; ii < K_factor; ++ii)
					{
						state_temp[0]  =  v_new[n];
						// state[(n * STATE_LEN) + 0] = v_new[n];
						// double Itot = TNNP_MarkovIKr_ODE(dt, - stims[n], &state[n * STATE_LEN], 1.0 , cell_type_vec.at(n));
						// double Itot = CNZ_ODE(dt, -stims, cell_para[n], &state[n * STATE_LEN]);
						// double Itot = INa_IKur_Drug_Updated_CNZ_ODE(dt, -stims, cell_para[n], &state[n * STATE_LEN]);
						double Itot = NCZ_ODE_Update(small_dt, -stims, cell_para[n], &state_temp[0]);

						v_new[n] = v_new[n] -  small_dt * Itot;
					}

					for (int jj = 1; jj < STATE_LEN; ++jj) {
						state[(n * STATE_LEN) + jj] = state_temp[jj];   // make a copy of state variables
					}
				} else {

					state[(n * STATE_LEN) + 0] = v_new[n];
					// double Itot = TNNP_MarkovIKr_ODE(dt, - stims[n], &state[n * STATE_LEN], 1.0 , cell_type_vec.at(n));
					// double Itot = CNZ_ODE(dt, -stims, cell_para[n], &state[n * STATE_LEN]);
					// double Itot = INa_IKur_Drug_Updated_CNZ_ODE(dt, -stims, cell_para[n], &state[n * STATE_LEN]);
					double Itot = NCZ_ODE_Update(dt, -stims, cell_para[n], &state[n * STATE_LEN]);

					if (fabs(Itot) <= 1)
						v_new[n] = v_new[n] -  dt * Itot;
					else {

						K_factor = 1 + (int) (fabs(Itot));
						K_factor = K_factor > 5 ? 5 : K_factor;
						double small_dt = dt / K_factor;
						for (int ii = 0; ii < K_factor; ++ii)
						{
							state_temp[0]  =  v_new[n];
							// state[(n * STATE_LEN) + 0] = v_new[n];
							// double Itot = TNNP_MarkovIKr_ODE(dt, - stims[n], &state[n * STATE_LEN], 1.0 , cell_type_vec.at(n));
							// double Itot = CNZ_ODE(dt, -stims, cell_para[n], &state[n * STATE_LEN]);
							// double Itot = INa_IKur_Drug_Updated_CNZ_ODE(dt, -stims, cell_para[n], &state[n * STATE_LEN]);
							double Itot = NCZ_ODE_Update(small_dt, -stims, cell_para[n], &state_temp[0]);

							v_new[n] = v_new[n] -  small_dt * Itot;
						}

						for (int jj = 1; jj < STATE_LEN; ++jj) {
							state[(n * STATE_LEN) + jj] = state_temp[jj];   // make a copy of state variables
						}
					}
				}


				/*if (fabs(Itot) < 1) {
					v_new[n] = v_new[n] -  dt * Itot;
				}*/

				// double Itot = TentVentByJohn(dt, -stims[n], &state[n * STATE_LEN], 1.0, cell_type_vec.at(n));
				// v_new[n] = cell_vec.at(n).svolt;

				// remember to update the voltage in the state array!!!! before proceed to the
				// state[(n * STATE_LEN) + 0] = v_new[n];
				// double Itot = crn_con_Vent(dt, &state[n * STATE_LEN], /*Zindex*/0, tissue[n]);
				// cell_vec.at(n).SetIStim();
				if (v_new[n] != v_new[n]) {
					std::cerr << " NaNs Encountered, Exiting Programe... ... \n\n\n";
					std::exit(0);
				}

				state[(n * STATE_LEN)] = v_new[n];
			}


#ifdef PROFILE_TIME
			#pragma omp single
			{
				// end_setup = omp_get_wtime();
				// fprintf(stderr, "Setup time:   %6g\n", end_setup - start);
				std::cerr << "Setup time:  " <<  omp_get_wtime() - end_setup << std::endl;
				timeelasped = omp_get_wtime();
			}
#endif

			/*operator spliting part 2*/

			/* using Heun's Method -> https://en.wikipedia.org/wiki/Heun%27s_method */
			/*#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_new[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				dv_diff_temp[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_new[n] + sum * dt / 2.0;  // devide by two because of operator splitting.
			}

			#pragma omp for schedule(static)   // part 2, y_i+1 = y_i + h/2*(f(t, y(t)) + f(t+1, y_bar_i+1))
			for (int n = 0; n < count; n++) {
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;
				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]);  // this is f(t+1, y_bar_i+1)
				}
				// dv[n] = sum;
				//y_bar_i+1
				v_old[n] = v_old[n] + (dt / 2.0) / 2.0 * (dv_diff_temp[n] + sum);  // devide by two because of operator splitting.
				APInfor_vec[n]->MeasureAPD90Using_dVdtMax(t, 0, Config_para.BCL, dt, v_old[n]);  // measure APD at the end of each time step.
			}*/




			/* using RK4 here*/
			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_new[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K1[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_new[n] + sum * dt / 2.0 / 2.0; // devide by two because of operator splitting. // y_n_K1/2
			}
			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K2[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp_2[n] = v_new[n] + sum * dt / 2.0 / 2.0; // devide by two because of operator splitting. // y_n_K2/2
			}

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp_2[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K3[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_new[n] + sum * dt / 2.0; // devide by two because of operator splitting. // y_n_K3/2
			}

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				// RK_K3[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_old[n] = v_old[n] + (sum + RK_K1[n] + 2 * RK_K2[n] + 2 * RK_K3[n]) * dt / 2.0 / 6.0; // devide by two because of operator splitting. // y_n_K3/2
				APInfor_vec[n]->MeasureAPD90Using_dVdtMax(t, 0, Config_para.BCL, dt, v_old[n]);  // measure APD at the end of each time step.
			}




			/* Euler method applied 4 times */
			/*#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_new[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K1[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_new[n] + sum * dt / 2.0 / 4.0; // devide by two because of operator splitting. // y_n_K1/2
			}
			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K2[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp_2[n] = v_temp[n] + sum * dt / 2.0 / 4.0; // devide by two because of operator splitting. // y_n_K2/2
			}

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp_2[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				RK_K3[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_temp[n] = v_temp_2[n] + sum * dt / 2.0/4.0; // devide by two because of operator splitting. // y_n_K3/2
			}

			#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				// v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_temp[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				// RK_K3[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_old[n] = v_temp[n] + sum * dt/2.0/4.0; // + (sum + RK_K1[n] + 2 * RK_K2[n] + 2 * RK_K3[n]) * dt / 2.0 / 6.0; // devide by two because of operator splitting. // y_n_K3/2
				APInfor_vec[n]->MeasureAPD90Using_dVdtMax(t, 0, Config_para.BCL, dt, v_old[n]);  // measure APD at the end of each time step.
			}
			*/


			/*Forward Euler */

			/*#pragma omp for schedule(static)   // part 1, y_bar_j+1 = y_i + h*f(t, y(t))
			for (int n = 0; n < count; n++) {

				v_old[n] = v_new[n];  // keep a copy of v_new[n], which is y_i;
				//dv[n] = 0.0;
				double sum = 0.0;
				int ii = 0;

				// #pragma omp simd reduction // not available on icpc 12.0
				for (ii = 0; ii < NUM_NBH; ii++) {
					sum +=  (lap[(n * NUM_NBH) + ii] * v_new[nbd[(n * NUM_NBH) + ii]]); // this is f(t, y_i)
				}
				// dv[n] = sum;
				//y_bar_i+1
				dv_diff_temp[n] = sum; // clear to 0 first, otherwise its going to accumulate ...
				v_old[n] = v_new[n] + sum * dt / 2.0;  // devide by two because of operator splitting.
				APInfor_vec[n]->MeasureAPD90Using_dVdtMax(t, 0, Config_para.BCL, dt, v_old[n]);  // measure APD at the end of each time step.
			}*/





#ifdef CV_COM
			#pragma omp for schedule(static)
			for (int n = 0; n < count; ++n) {

				/*#ifdef OUT_INICOND

								if (t >= Config_para.Total_time - 4*Config_para.BCL - 10 - dt / 2.0 and t <= Config_para.Total_time - 4*Config_para.BCL - 10 + dt / 2.0) // output initial values
								{
									char filename [100];
									sprintf(filename, "IniCond/Cell_%03d_BCL_%04d.bin", n, (int) Config_para.BCL);
									output_double_array_bin(filename, &state[n * STATE_LEN], STATE_LEN);
								}

				#endif
				*/

				PrePotential[n] = TemPotential[n];
				TemPotential[n] = v_old[n];
				// v_new[n] = v_new[n] -  dt * (Itot - stims[n]);
				if ((TemPotential[n] > -20.0) and (PrePotential[n] <= -20.0)) {

					if (t - ActTime[n] >= Config_para.S2 / 5.0)
					{
						ActTime[n] = t;
						AP_Number[n] += 1; // add An AP count
					}

				}




			}

			#pragma omp single
			{
				// printf("a\n");
				if ((TemPotential[FIRST_ACT] > -20.0) && (PrePotential[FIRST_ACT] <= -20.0)) {
					t_first_act =  t;
				}
				if ((TemPotential[SECOND_ACT] > -20.0) && (PrePotential[SECOND_ACT] <= -20.0)) {
					t_second_act =  t;
					double cvvalue = (SECOND_ACT - FIRST_ACT) * dx / (ActTime[SECOND_ACT] - ActTime[FIRST_ACT]);
					std::cerr << "t_first_act = " << t_first_act << " , t_second_act = " << t_second_act << std::endl;
					std::cerr << "CV = " << cvvalue << std::endl;
					cv.push_back(cvvalue);
#ifdef BREAK_IF_CVM

					if (cv.size() > 20)
					{
						if (fabs(cv[cv.size() - 1] - cv[cv.size() - 2]) < 0.0005)
						{
							// break;
							Config_para.Total_time = t;
						}
					}
#endif
				}

				/*std::ofstream out35("35.dat", std::ios::out | std::ios::app);
				out35 << t << " " << TemPotential[9] << std::endl;
				out35.close();*/
			}
#endif






#ifdef PROFILE_TIME
			#pragma omp single
			{
				// end_setup = omp_get_wtime();
				// fprintf(stderr, "Setup time:   %6g\n", end_setup - start);
				std::cerr << "Setup time:  " <<  omp_get_wtime() - end_setup << std::endl;
				timeelasped = omp_get_wtime();
			}
#endif
			timeint = t * intstim;
		} // end of time loop
#ifdef OUT_TO_BIN
		if (outindex != 0)
			if (omp_get_thread_num() < outindex) {
				output_voltage_array(count, v_out, "v%04d.bin", (filecount * omp_get_num_threads()) + omp_get_thread_num());
			}

#endif
		#pragma omp single
		{
			// for (int i = 0; i < count; ++i)
			{
				std::ofstream ConductionFile("conduction.log", std::ios::out);
				ConductionFile << AP_Number[0] << std::endl;
				ConductionFile << AP_Number[NX - 1] << std::endl;
				ConductionFile.close();
			}
		}
	} // end parallel loops

	delete [] Fcell;

	end = omp_get_wtime();
	std::cerr << "Elapsed Time : " << end - start << std::endl;
	std::cerr << "Time per ms : " << (end - end_setup) / (t_max - Config_para.t_start) << std::endl;

	/*for (int i = cv.size() - 4; i < cv.size(); ++i)
	{
		std::cout << Config_para.BCL <<  " "
		          << cv[i] <<  " "
		          << AP_Number[0] << " "
		          << AP_Number[NX - 1] << std::endl;
	}*/
	// if (cv.size() == 2)

	try {
		std::cout << Config_para.S2 << " " << cv.back() << " " << APInfor_vec[FIRST_ACT]->APD_Vec.back() << " " << APInfor_vec[SECOND_ACT]->APD_Vec.back() << " "

		          << APInfor_vec[SECOND_ACT]->dVdtmax_Vec.back()  << " "
		          <<  AP_Number[0] << " "
		          << AP_Number[NX - 1] << " "
		          << Config_para.S2 << " " << cv[cv.size() - 2] << " " << APInfor_vec[FIRST_ACT]->APD_Vec[APInfor_vec[FIRST_ACT]->APD_Vec.size() - 2] << " "
		          << APInfor_vec[SECOND_ACT]->APD_Vec[APInfor_vec[SECOND_ACT]->APD_Vec.size() - 2] << " "
		          << APInfor_vec[SECOND_ACT]->dVdtmax_Vec [APInfor_vec[SECOND_ACT]->dVdtmax_Vec.size() - 2] << " "
		          <<  AP_Number[0] << " "
		          << AP_Number[NX - 1] << " " << std::endl;
	} catch (std::exception & e) {
		std::cerr << "An exception occurred. Exception Nr. " << e.what() << " at the very end of the programme " << '\n';
	}


	// if want to output to see all APDs in each cell.
#ifdef OUT_ALL_APD
	try {

		for (int i = 0; i < APInfor_vec[0]->APD_Vec.size(); ++i)
		{
			for (auto& x : APInfor_vec)
				std::cout << x->APD_Vec[i] << " " ;
			std::cout << std::endl;
		}

	} catch (std::exception & e) {
		std::cerr << "An exception occurred. Exception Nr. " << e.what() << " at the very end of the programme " << '\n';
	}

#endif




	std::ofstream OneD_APD_out;
	OneD_APD_out.open(Config_para.OneD_OutFile, std::ios::app);


	int n1 = FIRST_ACT;
	// int end = cv.size() - 1;
	// int end_2 = cv.size() - 2;
	int n2 = SECOND_ACT;

	OneD_APD_out << Config_para.Sim_ID << " "
	             << cv.back() << " "
	             << cv[cv.size() - 2] << " "
	             << AP_Number[0] << " "
	             << AP_Number[NX - 1] << " "
	             << Config_para.BCL << " "
	             << APInfor_vec[n1]->APD_out_end[0] << " "
	             << APInfor_vec[n1]->APD75_out[0] << " "
	             << APInfor_vec[n1]->APD50_out[0] << " "
	             << APInfor_vec[n1]->APD30_out[0] << " "
	             << APInfor_vec[n1]->APD20_out[0] << " "
	             << APInfor_vec[n1]->Vmax_out[0] << " "
	             << APInfor_vec[n1]->Vmin_out[0] << " "
	             << APInfor_vec[n1]->dVdtmax_out[0] << " "
	             << APInfor_vec[n1]->plateau_potential_out[0] << " "
	             << APInfor_vec[n1]->INa_max_out[0] << " "
	             << APInfor_vec[n1]->dVdt_Time_RP[0] << " "
	             << APInfor_vec[n1]->V_m70_Time_RP[0] << " "
	             << APInfor_vec[n1]->Cai_max_out[0] << " "
	             << APInfor_vec[n1]->Cai_min_out[0] << " "
	             /* std::cout*/ << Config_para.BCL << " "
	             << APInfor_vec[n1]-> APD_out_end[1] << " "
	             << APInfor_vec[n1]-> APD75_out[1] << " "
	             << APInfor_vec[n1]-> APD50_out[1] << " "
	             << APInfor_vec[n1]-> APD30_out[1] << " "
	             << APInfor_vec[n1]-> APD20_out[1] << " "
	             << APInfor_vec[n1]-> Vmax_out[1] << " "
	             << APInfor_vec[n1]-> Vmin_out[1] << " "
	             << APInfor_vec[n1]-> dVdtmax_out[1] << " "
	             << APInfor_vec[n1]-> plateau_potential_out[1] << " "
	             << APInfor_vec[n1]-> INa_max_out[1] << " "
	             << APInfor_vec[n1]-> dVdt_Time_RP[1] << " "
	             << APInfor_vec[n1]-> V_m70_Time_RP[1] << " "
	             << APInfor_vec[n1]-> Cai_max_out[1] << " "
	             << APInfor_vec[n1]-> Cai_min_out[1] << " "
	             << Config_para.BCL << " "
	             << APInfor_vec[n2]->APD_out_end[0] << " "
	             << APInfor_vec[n2]->APD75_out[0] << " "
	             << APInfor_vec[n2]->APD50_out[0] << " "
	             << APInfor_vec[n2]->APD30_out[0] << " "
	             << APInfor_vec[n2]->APD20_out[0] << " "
	             << APInfor_vec[n2]->Vmax_out[0] << " "
	             << APInfor_vec[n2]->Vmin_out[0] << " "
	             << APInfor_vec[n2]->dVdtmax_out[0] << " "
	             << APInfor_vec[n2]->plateau_potential_out[0] << " "
	             << APInfor_vec[n2]->INa_max_out[0] << " "
	             << APInfor_vec[n2]->dVdt_Time_RP[0] << " "
	             << APInfor_vec[n2]->V_m70_Time_RP[0] << " "
	             << APInfor_vec[n2]->Cai_max_out[0] << " "
	             << APInfor_vec[n2]->Cai_min_out[0] << " "
	             /* std::cout*/ << Config_para.BCL << " "
	             << APInfor_vec[n2]-> APD_out_end[1] << " "
	             << APInfor_vec[n2]-> APD75_out[1] << " "
	             << APInfor_vec[n2]-> APD50_out[1] << " "
	             << APInfor_vec[n2]-> APD30_out[1] << " "
	             << APInfor_vec[n2]-> APD20_out[1] << " "
	             << APInfor_vec[n2]-> Vmax_out[1] << " "
	             << APInfor_vec[n2]-> Vmin_out[1] << " "
	             << APInfor_vec[n2]-> dVdtmax_out[1] << " "
	             << APInfor_vec[n2]-> plateau_potential_out[1] << " "
	             << APInfor_vec[n2]-> INa_max_out[1] << " "
	             << APInfor_vec[n2]-> dVdt_Time_RP[1] << " "
	             << APInfor_vec[n2]-> V_m70_Time_RP[1] << " "
	             << APInfor_vec[n2]-> Cai_max_out[1] << " "
	             << APInfor_vec[n2]-> Cai_min_out[1] << " "

	             << std::endl;

	OneD_APD_out.close();


	return 1;
}
/* end of main() */
