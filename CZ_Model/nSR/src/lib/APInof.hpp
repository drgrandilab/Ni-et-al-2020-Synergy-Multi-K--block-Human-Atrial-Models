/*

    Update: accurate measurement of RMP. (or MDP most distolic Potential)

    Haibo Ni
    Fri 26 Aug 2016 17:51:50 BST


    updated: calculate plateau potentials by Haibo Ni Wed 06 Jul 2016 15:45:53 BST

    simple function to measure the APDs

    Haibo Ni
    qiangzi.ni@gmail.com
    update Jan. 27. 2015
    Wed 25 May 2016 15:10:59 BST
*/


#ifndef AP_INFOR_H
#define AP_INFOR_H

#include <iostream>
#include <fstream>
#include <cmath>


class APInfor {
public:
    double Vmax, Vmin, timeAPDstart, timeAPDend, dVdtmax, APD_out_end[2], Vmin_last;
    double APD50_out[2], APD30_out[2], APD75_out[2], APD20_out[2];
    double Vmax_out[2], Vmin_out[2], dVdtmax_out[2], INa_max_out[2];
    double plateau_potential_out[2];
    int APD_switch, APD_count, Vmax_switch, APD_out_swhich;
    double APD90, APD50, APD30, APD20, APD75;
    double Vm_prev;
    double AMP, AMP_last, AMP_over_last;
    double APD90_prev;
    double t_since_up_stroke;
    double v90, v75, v50, v30, v20;
    double ICaL_in, RyR_in;
    double Istim_prev;
    double Diastolic_value;
    double Systolic_value_H, Systolic_value_L;
    unsigned int N_stim;
    double Diastolic_value_2, Systolic_value_L_2, Systolic_value_H_2;
    double INa_max;
    double plateau_potential;
    double dVdt_Time, V_m70_Time;
    double dVdt_Time_RP[2];
    double V_m70_Time_RP[2];
    double Cai_max, Cai_min, Cai_max_out[2], Cai_min_out[2];


    std::ofstream out;

    APInfor(const char *filename  = "APDmeasure.dat", bool file_app_mode = false);
    ~APInfor();

    void MeasureAPD90(double t, double Istim, double BCL, double dt, double Vm);
    void MeasureAPD90_INa(double t, double Istim, double BCL, double dt, double Vm, double INa, double value);
    void MeasureAPD90andDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value);
    void MeasureAPD90andTwoDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value_1, double Value_2);
    void MeasureAPD90andDSValuewith_INa(double INa_threshold, double INa, double t, double Istim, double BCL, double dt, double Vm, double Value);
    void MeasureAPD90andDSValuewithStrokeTime(double StrokeTime, double t, double Istim, double BCL, double dt, double Vm, double Value);
    void CalciumAccumulation(double t, double Istim, double BCL, double dt, double ICaL_con, double J_RyR);
    void ReportAPD();
    void ReportLastTwo();
    void ReportLastTwo(double);
    void ReportLast();
};


APInfor::APInfor(const char *filename, bool file_app_mode)  {
    Vmax = 0.0;
    Vmin = 0.0;
    timeAPDstart = 0.0;
    timeAPDend = 0.0;
    dVdtmax = 0.0;
    APD_out_end[0] = 0.0;
    APD_out_end[1] = 0.0;

    APD_switch = 0;
    APD_count = 0;
    Vmax_switch = 0;
    APD_out_swhich = 0;
    APD90 = APD75 = APD50 = APD30 = APD20 = 0.0;
    Vm_prev = 0.0;
    v90 = 0.0;
    v75 = 0.0;
    v50 = 0.0;
    v30 = 0.0;
    v20 = 0.0;
    ICaL_in = 0.0;
    RyR_in = 0.0;
    Istim_prev = 0.0;
    AMP = 0.0;
    AMP_last = 0.0;;
    AMP_over_last = 0.0;
    APD90_prev = 0.0;
    t_since_up_stroke = 0.0;
    plateau_potential = 0.0;
    Diastolic_value = 0.0;
    Systolic_value_L = 0.0;
    Systolic_value_H = 0.0;
    Diastolic_value_2 = 0.0;
    Systolic_value_L_2 = 0.0;
    Systolic_value_H_2 = 0.0;
    N_stim = 0;
    INa_max = 0.0;
    Vmin_last = 0.0;
    dVdt_Time = 0;
    /* are you going to output in appending mode? */
    if (file_app_mode)
    {
        out.open(filename, std::ios::out | std::ios::app);
    } else {
        out.open(filename, std::ios::out);
    }
}

APInfor::~APInfor() {
    if (out.is_open()) {
        out.close();
    }
}
void APInfor::MeasureAPD90(double t, double Istim, double BCL, double dt, double Vm)  {


    if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
        APD_switch = 1;
        Vmax = Vm;
        timeAPDstart = t;
        // Vmin = Vm;

        if (Vmin > Vm)
        {
            Vmin = Vm; // for the first Beat... added Haibo
        }
        dVdtmax = 0.0;
        N_stim ++;
        Vmax_switch = 0;
        INa_max = 0.0;
        plateau_potential = 0.0;
    }


    if (APD_switch == 0)
    {
        if (Vmin >= Vm) {
            Vmin = Vm;
        }
    }


    if (APD_switch == 1) {
        if (Vm > -60.0) {
            if (Vm >= Vmax) {
                Vmax = Vm;
            }
            else Vmax_switch = 1;
            if (Vmax_switch == 1) {
                APD_switch = 2;
                AMP_last = AMP;
                AMP = Vmax - Vmin;
                AMP_over_last = AMP / AMP_last;
                /*printf("%f %f\n", Vmax, Vmin);
                printf("%f %f\n", AMP, AMP_last);*/
                v90 = Vmax - 0.9 * (Vmax - Vmin);
                v75 = Vmax - 0.75 * (Vmax - Vmin);
                v50 = Vmax - 0.50 * (Vmax - Vmin);
                v30 = Vmax - 0.30 * (Vmax - Vmin);

                Vmin_last = Vmin;
                Vmin = 0;
                // printf("asfde\n");
            }
        }
    }

    /*

        plateau_potential is defined as the average potential
       between 10 and 50 ms following application of stimulation.
       Haibo Ni.
    */

    if (t > timeAPDstart + 10 and t < timeAPDstart + 50 )
    {
        plateau_potential = plateau_potential + dt * Vm;
    }

    if (APD_switch == 2) {

        if ((Vm_prev >= v30) && (Vm <= v30) ) {
            APD30 = t - timeAPDstart ;
        }
        else if ((Vm_prev >= v50) && (Vm <= v50 )) {
            APD50 = t - timeAPDstart ;
        }
        else if (Vm_prev >= v75 && Vm <= v75 ) {
            APD75 = t - timeAPDstart ;
        }
        else if (Vm_prev >= v90 && Vm <= v90) {
            APD_switch = 0;
            APD_count ++;
            Vmax_switch = 0;
            APD90 = t - timeAPDstart;
            plateau_potential /= (50.0 - 10.0); // average
            // printf("dddddddddddd\n");
            out << APD_count *BCL << " "
                << APD90 << " "
                << APD75 << " "
                << APD50 << " "
                << APD30 << " "
                << Vmax  << " "
                << Vmin_last  << " "
                << dVdtmax  << " "
                << plateau_potential  << " "
                << AMP_over_last  << " "
                << std::endl;

            if (APD_out_swhich == 0) {
                APD_out_end[0] = t - timeAPDstart;
                APD75_out[0] = APD75;
                APD50_out[0] = APD50;
                APD30_out[0] = APD30;
                Vmax_out[0] = Vmax;
                Vmin_out[0] = Vmin_last;
                dVdtmax_out[0] = dVdtmax;
                APD_out_swhich = 1;
                plateau_potential_out[0] = plateau_potential;
            }
            else if (APD_out_swhich == 1) {
                APD_out_end[1] = t - timeAPDstart;
                APD75_out[1]   = APD75;
                APD50_out[1]   = APD50;
                APD30_out[1]   = APD30;
                Vmax_out[1]    = Vmax;
                Vmin_out[1]    = Vmin_last;
                dVdtmax_out[1] = dVdtmax;
                plateau_potential_out[1] = plateau_potential;
                APD_out_swhich = 0;
            }
        }
    }

    double dVdt = (Vm - Vm_prev) / dt;
    if (dVdt > dVdtmax) dVdtmax = dVdt;

    Vm_prev = Vm;
    Istim_prev = Istim;
}


void APInfor::MeasureAPD90_INa(double t, double Istim, double BCL, double dt, double Vm, double INa, double value) {

    if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {


        out << std::setprecision(9) << APD_count *BCL << " "
            << APD90 << " "
            << APD75 << " "
            << APD50 << " "
            << APD30 << " "
            << APD20 << " "
            << Vmax  << " "
            << Vmin_last  << " "
            << dVdtmax  << " "
            << plateau_potential  << " "
            << INa_max << " "
            << dVdt_Time  << " "
            << V_m70_Time << " "
            << AMP_over_last  << " "
            << std::endl;
        if (APD_out_swhich == 0) {
            APD_out_end[0] = APD90;
            APD75_out[0] = APD75;
            APD50_out[0] = APD50;
            APD30_out[0] = APD30;
            APD20_out[0] = APD20;
            Vmax_out[0] = Vmax;
            Vmin_out[0] = Vmin_last;
            dVdtmax_out[0] = dVdtmax;
            plateau_potential_out[0] = plateau_potential;
            INa_max_out[0] = INa_max;
            dVdt_Time_RP[0] = dVdt_Time;
            V_m70_Time_RP[0] = V_m70_Time;
            Cai_max_out[0] = Cai_max;
            Cai_min_out[0] = Cai_min;
            APD_out_swhich = 1;
        }
        else if (APD_out_swhich == 1) {
            APD_out_end[1] = APD90;
            APD75_out[1]   = APD75;
            APD50_out[1]   = APD50;
            APD30_out[1]   = APD30;
            APD20_out[1]   = APD20;
            Vmax_out[1]    = Vmax;
            Vmin_out[1]    = Vmin_last;
            dVdtmax_out[1] = dVdtmax;
            INa_max_out[1] = INa_max;
            dVdt_Time_RP[1] = dVdt_Time;
            V_m70_Time_RP[1] = V_m70_Time;
            plateau_potential_out[1] = plateau_potential;
            Cai_max_out[1] = Cai_max;
            Cai_min_out[1] = Cai_min;
            APD_out_swhich = 0;
        }


        APD_switch = 1;
        Vmax = Vm;
        timeAPDstart = t;
        // Vmin = Vm;

        if (Vmin > Vm)
        {
            Vmin = Vm; // for the first Beat... added Haibo
        }

        Cai_max = -100;
        Cai_min = 100;
        dVdtmax = 0.0;
        N_stim ++;
        Vmax_switch = 0;
        INa_max = 0.0;
        plateau_potential = -100.0;
        APD90 = -1.0;
        APD75 = -1.0;
        APD50 = -1.0;
        APD30 = -.0;
        APD20 = -1.0;
        Vmax = -100;
        Vmin_last = 0.0;
        dVdtmax = 0.0;
        INa_max = 0.0;
        dVdt_Time = -1;
        V_m70_Time = -1;
    }


    if (APD_switch == 0)
    {
        if (Vmin >= Vm) {
            Vmin = Vm;
        }
    }



    if (Cai_max < value)
    {
        Cai_max = value;
    }
    if (Cai_min > value)
    {
        Cai_min = value;
    }


    if (APD_switch == 2) {

        if (Vm >= Vmax) {
            Vmax = Vm;
            // APD_switch =1;
        }
    }


    if (APD_switch == 1) {

        if (INa_max > INa)
        {
            INa_max = INa;
            // std::cout << "yes" << std::endl;
        }

        if (Vm > -75.0) {

            if (Vm >= Vmax) {
                Vmax = Vm;
            }
            else Vmax_switch = 1;
            if (Vmax_switch == 1) {
                APD_switch = 2;
                AMP_last = AMP;
                AMP = Vmax - Vmin;
                AMP_over_last = AMP / AMP_last;
                /*printf("%f %f\n", Vmax, Vmin);
                printf("%f %f\n", AMP, AMP_last);*/
                v90 = Vmax - 0.9 * (Vmax - Vmin);
                v75 = Vmax - 0.75 * (Vmax - Vmin);
                v50 = Vmax - 0.50 * (Vmax - Vmin);
                v30 = Vmax - 0.30 * (Vmax - Vmin);
                v20 = Vmax - 0.20 * (Vmax - Vmin);
                Vmin_last = Vmin;
                Vmin = 0;
                // printf("asfde\n");
            }
        }
    }

    /*

        plateau_potential is defined as the average potential
       between 10 and 50 ms following application of stimulation.
       Haibo Ni.
    */

    if (t > timeAPDstart + 10 and t < timeAPDstart + 50 )
    {
        plateau_potential = plateau_potential + dt * Vm;
    }

    if (APD_switch == 2) {

        if (fabs(Vm - (-70)) < 0.1) {
            V_m70_Time = t;
        }
        if ((Vm_prev >= v20) && (Vm <= 20) ) {
            APD20 = t - timeAPDstart ;
        } else if ((Vm_prev >= v30) && (Vm <= v30) ) {
            APD30 = t - timeAPDstart ;
        }
        else if ((Vm_prev >= v50) && (Vm <= v50 )) {
            APD50 = t - timeAPDstart ;
        }
        else if (Vm_prev >= v75 && Vm <= v75 ) {
            APD75 = t - timeAPDstart ;
        }
        else if (Vm_prev >= v90 && Vm <= v90) {
            APD_switch = 0;
            APD_count ++;
            Vmax_switch = 0;
            APD90 = t - timeAPDstart;
            plateau_potential /= (50.0 - 10.0); // average
            // printf("dddddddddddd\n");

        }
    }
    double dVdt = (Vm - Vm_prev) / dt;
    if (dVdt > dVdtmax) {
        dVdtmax = dVdt;
        dVdt_Time = t;
    }

    Vm_prev = Vm;
    Istim_prev = Istim;
}



// void APInfor::MeasureAPD90_INa(double t, double Istim, double BCL, double dt, double Vm, double INa)  {

//     if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
//         APD_switch = 1;
//         Vmax = Vm;
//         timeAPDstart = t;
//         Vmin = Vm;
//         dVdtmax = 0.0;
//         N_stim ++;
//         Vmax_switch = 0;
//         INa_max = 0.0;
//     }

//     if (APD_switch == 1) {
//         if (Vm > -30.0) {

//             if (INa_max > INa)
//             {
//                 INa_max = INa;
//                 // std::cout << "yes" << std::endl;
//             }


//             if (Vm >= Vmax) {
//                 Vmax = Vm;
//             }
//             else Vmax_switch = 1;
//             if (Vmax_switch == 1) {
//                 APD_switch = 2;
//                 AMP_last = AMP;
//                 AMP = Vmax - Vmin;
//                 AMP_over_last = AMP / AMP_last;
//                 /*printf("%f %f\n", Vmax, Vmin);
//                 printf("%f %f\n", AMP, AMP_last);*/
//                 v90 = Vmax - 0.9 * (Vmax - Vmin);
//                 v75 = Vmax - 0.75 * (Vmax - Vmin);
//                 v50 = Vmax - 0.50 * (Vmax - Vmin);
//                 v30 = Vmax - 0.30 * (Vmax - Vmin);
//                 // printf("asfde\n");
//             }
//         }
//     }
//     if (APD_switch == 2) {

//         if ((Vm_prev >= v30) && (Vm <= v30) ) {
//             APD30 = t - timeAPDstart ;
//         }
//         else if ((Vm_prev >= v50) && (Vm <= v50 )) {
//             APD50 = t - timeAPDstart ;
//         }
//         else if (Vm_prev >= v75 && Vm <= v75 ) {
//             APD75 = t - timeAPDstart ;
//         }
//         else if (Vm_prev >= v90 && Vm <= v90) {
//             APD_switch = 0;
//             APD_count ++;
//             Vmax_switch = 0;
//             APD90 = t - timeAPDstart;
//             // printf("dddddddddddd\n");
//             if (t > 20000) {

//                 out << APD_count *BCL << " "
//                     << APD90 << " "
//                     << APD75 << " "
//                     << APD50 << " "
//                     << APD30 << " "
//                     << Vmax  << " "
//                     << Vmin  << " "
//                     << dVdtmax << " "
//                     << INa_max << " "
//                     << BCL << " "
//                     << std::endl;
//             }

//             if (APD_out_swhich == 0) {
//                 APD_out_end[0] = t - timeAPDstart;
//                 APD75_out[0] = APD75;
//                 APD50_out[0] = APD50;
//                 APD30_out[0] = APD30;
//                 Vmax_out[0] = Vmax;
//                 Vmin_out[0] = Vmin;
//                 dVdtmax_out[0] = dVdtmax;
//                 APD_out_swhich = 1;
//             }
//             else if (APD_out_swhich == 1) {
//                 APD_out_end[1] = t - timeAPDstart;
//                 APD75_out[1]   = APD75;
//                 APD50_out[1]   = APD50;
//                 APD30_out[1]   = APD30;
//                 Vmax_out[1]    = Vmax;
//                 Vmin_out[1]    = Vmin;
//                 dVdtmax_out[1] = dVdtmax;
//                 APD_out_swhich = 0;
//             }
//         }
//     }

//     double dVdt = (Vm - Vm_prev) / dt;
//     if (dVdt > dVdtmax) dVdtmax = dVdt;

//     Vm_prev = Vm;
//     Istim_prev = Istim;
// }




void APInfor::MeasureAPD90andDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value)  {

    if (((Istim >= 1e-8 || Istim <= -1e-8 ) && (fabs(Istim_prev) < 1e-10) && (APD_switch == 0))) {
        if (APD_count > 0)
        {

            out << APD_count *BCL << " "
                << BCL << " "
                << APD90 << " "
                << APD75 << " "
                << APD50 << " "
                << APD30 << " "
                << Systolic_value_L << " "
                << Systolic_value_H << " "

                << Vmax << " "
                << Vmin << " "
                << dVdtmax << std::endl;

        }
        APD_switch = 1;
        Vmax = -60;
        timeAPDstart = t;
        Vmin = Vm;
        dVdtmax = 0.0;
        Diastolic_value = Value;
        Systolic_value_L = Value;
        Systolic_value_H = Value;
        N_stim ++;
        Vmax_switch = 0;

        // printf("ssssssssss\n");
    }
    if (Systolic_value_L > Value)
    {
        Systolic_value_L = Value;
    }
    if (Systolic_value_H < Value)
    {
        Systolic_value_H = Value;
    }
    if (APD_switch == 1) {
        if (Vm > -80.0) {
            if (Vm >= Vmax) {
                Vmax = Vm;
            }
            else Vmax_switch = 1;
            if (Vmax_switch == 1) {
                APD_switch = 2;
                AMP_last = AMP;
                AMP = Vmax - Vmin;
                AMP_over_last = AMP / AMP_last;
                /*printf("%f %f\n", Vmax, Vmin);

                printf("%f %f\n", AMP, AMP_last);*/
                v90 = Vmax - 0.9 * (Vmax - Vmin);
                v75 = Vmax - 0.75 * (Vmax - Vmin);
                v50 = Vmax - 0.50 * (Vmax - Vmin);
                v30 = Vmax - 0.30 * (Vmax - Vmin);
                // printf("asfde\n");
            }
        }
    }
    if (APD_switch == 2) {

        if ((Vm_prev >= v30) && (Vm <= v30) ) {
            APD30 = t - timeAPDstart ;
            // printf("aaaa\n");

        }
        else if ((Vm_prev >= v50) && (Vm <= v50 )) {
            APD50 = t - timeAPDstart ;
            // printf("bbbbbbbbb\n");

        }
        else if (Vm_prev >= v75 && Vm <= v75 ) {
            APD75 = t - timeAPDstart ;
            // printf("ccccccccccc\n");

        }
        else if (Vm_prev >= v90 && Vm <= v90) {
            APD_switch = 0;
            APD_count ++;
            Vmax_switch = 0;
            APD90 = t - timeAPDstart;
            // printf("dddddddddddd\n");
            // printf("%f %f %f\n", Diastolic_value, Systolic_value_L, Systolic_value_H);

            // fprintf (out, "%.2f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", APD_count * BCL, APD90, APD75, APD50, APD30, Vmax, Vmin, dVdtmax);
            if (APD_out_swhich == 0) {
                APD_out_end[0] = t - timeAPDstart;
                APD_out_swhich = 1;
            }
            else if (APD_out_swhich == 1) {
                APD_out_end[1] = t - timeAPDstart;
                APD_out_swhich = 0;
            }
        }
    }

    double dVdt = (Vm - Vm_prev) / dt;
    if (dVdt > dVdtmax) dVdtmax = dVdt;
    Istim_prev = Istim;
    Vm_prev = Vm;
}
void APInfor::MeasureAPD90andTwoDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value_1, double Value_2)  {

    if (((Istim >= 1e-8 || Istim <= -1e-8 ) && (fabs(Istim_prev) < 1e-10) && (APD_switch == 0))) {
        if (APD_count > 70)
        {
            out << APD_count *BCL << " "
                << BCL << " "
                << APD90 << " "
                << APD75 << " "
                << APD50 << " "
                << APD30 << " "
                << Systolic_value_L << " "
                << Systolic_value_H << " "
                << Systolic_value_L_2 << " "
                << Systolic_value_H_2 << " "
                << Vmax << " "
                << Vmin << " "
                << dVdtmax << std::endl;
        }

        APD_switch = 1;
        Vmax = -60;
        timeAPDstart = t;
        Vmin = Vm;
        dVdtmax = 0.0;
        Diastolic_value = Value_1;
        Systolic_value_L = Value_1;
        Systolic_value_H = Value_1;
        Diastolic_value_2 = Value_2;
        Systolic_value_L_2 = Value_2;
        Systolic_value_H_2 = Value_2;
        N_stim ++;
        Vmax_switch = 0;

        // printf("ssssssssss\n");
    }
    if (Systolic_value_L > Value_1)
    {
        Systolic_value_L = Value_1;
    }
    if (Systolic_value_H < Value_1)
    {
        Systolic_value_H = Value_1;
    }
    if (Systolic_value_L_2 > Value_2)
    {
        Systolic_value_L_2 = Value_2;
    }
    if (Systolic_value_H_2 < Value_2)
    {
        Systolic_value_H_2 = Value_2;
    }
    if (APD_switch == 1) {
        if (Vm > -30.0) {
            if (Vm >= Vmax) {
                Vmax = Vm;
            }
            else Vmax_switch = 1;
            if (Vmax_switch == 1) {
                APD_switch = 2;
                AMP_last = AMP;
                AMP = Vmax - Vmin;
                AMP_over_last = AMP / AMP_last;
                /*printf("%f %f\n", Vmax, Vmin);

                printf("%f %f\n", AMP, AMP_last);*/
                v90 = Vmax - 0.9 * (Vmax - Vmin);
                v75 = Vmax - 0.75 * (Vmax - Vmin);
                v50 = Vmax - 0.50 * (Vmax - Vmin);
                v30 = Vmax - 0.30 * (Vmax - Vmin);
                // printf("asfde\n");
            }
        }
    }
    if (APD_switch == 2) {

        if ((Vm_prev >= v30) && (Vm <= v30) ) {
            APD30 = t - timeAPDstart ;
            // printf("aaaa\n");

        }
        else if ((Vm_prev >= v50) && (Vm <= v50 )) {
            APD50 = t - timeAPDstart ;
            // printf("bbbbbbbbb\n");

        }
        else if (Vm_prev >= v75 && Vm <= v75 ) {
            APD75 = t - timeAPDstart ;
            // printf("ccccccccccc\n");

        }
        else if (Vm_prev >= v90 && Vm <= v90) {
            APD_switch = 0;
            APD_count ++;
            Vmax_switch = 0;
            APD90 = t - timeAPDstart;
            // printf("dddddddddddd\n");
            // printf("%f %f %f\n", Diastolic_value, Systolic_value_L, Systolic_value_H);
            // fprintf (out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", APD_count * BCL, APD90, APD75, APD50, APD30, Vmax, Vmin, dVdtmax);
            if (APD_out_swhich == 0) {
                APD_out_end[0] = t - timeAPDstart;
                APD_out_swhich = 1;
            }
            else if (APD_out_swhich == 1) {
                APD_out_end[1] = t - timeAPDstart;
                APD_out_swhich = 0;
            }
        }
    }

    double dVdt = (Vm - Vm_prev) / dt;
    if (dVdt > dVdtmax) dVdtmax = dVdt;
    Istim_prev = Istim;
    Vm_prev = Vm;
}


/*INa is used to measure the upstroke */
void APInfor::MeasureAPD90andDSValuewith_INa(double INa_threshold, double INa, double t, double Istim, double BCL, double dt, double Vm, double Value)  {

    if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
        if (APD_count > 0)
        {
            // printf("dddddddddddd\n");
            // printf("%f %f %f\n", Diastolic_value, Systolic_value_L, Systolic_value_H);

            // fprintf (out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f  %.3f %.3f %.3f %.3f\n", APD_count * BCL, BCL, APD90, APD75, APD50, APD30, Systolic_value_L, Systolic_value_H, Vmax, Vmin, dVdtmax);
            out << APD_count *BCL << " "
                << BCL << " "
                << APD90 << " "
                << APD75 << " "
                << APD50 << " "
                << APD30 << " "
                << Systolic_value_L << " "
                << Systolic_value_H << " "
                << Vmax << " "
                << Vmin << " "
                << dVdtmax << std::endl;
        }
        APD_switch = 1;
        Vmax = -60;
        timeAPDstart = t;
        Vmin = Vm;
        dVdtmax = 0.0;
        Diastolic_value = Value;
        Systolic_value_L = Value;
        Systolic_value_H = Value;
        N_stim ++;
        Vmax_switch = 0;

        // printf("ssssssssss\n");
    }
    if (Systolic_value_L > Value)
    {
        Systolic_value_L = Value;
    }
    if (Systolic_value_H < Value)
    {
        Systolic_value_H = Value;
    }
    if (APD_switch == 1) {
        if (Vm > -30.0) {
            if ((Vm > Vmax) && (INa <= INa_threshold)) {
                Vmax = Vm;
            }
            else Vmax_switch = 1;
            if (Vmax_switch == 1) {
                APD_switch = 2;
                AMP_last = AMP;
                AMP = Vmax - Vmin;
                AMP_over_last = AMP / AMP_last;
                /*printf("%f %f\n", AMP, AMP_last);
                printf("%f %f\n", Vmax, Vmin);*/

                v90 = Vmax - 0.9 * (Vmax - Vmin);
                v75 = Vmax - 0.75 * (Vmax - Vmin);
                v50 = Vmax - 0.50 * (Vmax - Vmin);
                v30 = Vmax - 0.30 * (Vmax - Vmin);
                // printf("asfde\n");
            }
        }
    }
    if (APD_switch == 2) {

        if ((Vm_prev >= v30) && (Vm <= v30) ) {
            APD30 = t - timeAPDstart ;
            // printf("aaaa\n");

        }
        else if ((Vm_prev >= v50) && (Vm <= v50 )) {
            APD50 = t - timeAPDstart ;
            // printf("bbbbbbbbb\n");

        }
        else if (Vm_prev >= v75 && Vm <= v75 ) {
            APD75 = t - timeAPDstart ;
            // printf("ccccccccccc\n");

        }
        else if (Vm_prev >= v90 && Vm <= v90) {
            APD_switch = 0;
            APD_count ++;
            Vmax_switch = 0;
            APD90 = t - timeAPDstart;
            if (APD_out_swhich == 0) {
                APD_out_end[0] = t - timeAPDstart;
                APD_out_swhich = 1;
            }
            else if (APD_out_swhich == 1) {
                APD_out_end[1] = t - timeAPDstart;
                APD_out_swhich = 0;
            }
        }
    }

    double dVdt = (Vm - Vm_prev) / dt;
    if (dVdt > dVdtmax) dVdtmax = dVdt;
    Istim_prev = Istim;
    Vm_prev = Vm;
}

void APInfor::MeasureAPD90andDSValuewithStrokeTime(double StrokeTime, double t, double Istim, double BCL, double dt, double Vm, double Value)  {

    if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
        APD_switch = 1;
        Vmax = -60;
        timeAPDstart = t;
        Vmin = Vm;
        dVdtmax = 0.0;
        Diastolic_value = Value;
        Systolic_value_L = Value;
        Systolic_value_H = Value;
        N_stim ++;
        t_since_up_stroke = 0.0;
        Vmax_switch = 0;

        // printf("ssssssssss\n");
    }
    if (Systolic_value_L > Value)
    {
        Systolic_value_L = Value;
    }
    if (Systolic_value_H < Value)
    {
        Systolic_value_H = Value;
    }
    t_since_up_stroke += dt;

    if (APD_switch == 1) {

        if (Vm > -10.0) {
            if ((Vm > Vmax) && (t_since_up_stroke <= StrokeTime)) {
                Vmax = Vm;
                // printf("%f\n", t_since_up_stroke);
                // printf("%f\n", Vmax);

            }
            else Vmax_switch = 1;
            if (Vmax_switch == 1) {
                APD_switch = 2;
                AMP_last = AMP;
                AMP = Vmax - Vmin;
                AMP_over_last = AMP / AMP_last;
                // printf("%f %f\n", AMP, AMP_last);
                // printf("%f %f\n", Vmax, Vmin);

                v90 = Vmax - 0.9 * (Vmax - Vmin);
                v75 = Vmax - 0.75 * (Vmax - Vmin);
                v50 = Vmax - 0.50 * (Vmax - Vmin);
                v30 = Vmax - 0.30 * (Vmax - Vmin);
                // printf("asfde\n");
            }
        }
    }
    if (APD_switch == 2) {

        if ((Vm_prev >= v30) && (Vm <= v30) ) {
            APD30 = t - timeAPDstart ;
            // printf("aaaa\n");

        }
        else if ((Vm_prev >= v50) && (Vm <= v50 )) {
            APD50 = t - timeAPDstart ;
            // printf("bbbbbbbbb\n");

        }
        else if (Vm_prev >= v75 && Vm <= v75 ) {
            APD75 = t - timeAPDstart ;
            // printf("ccccccccccc\n");

        }
        else if (Vm_prev >= v90 && Vm <= v90) {
            APD_switch = 0;
            APD_count ++;
            Vmax_switch = 0;
            APD90 = t - timeAPDstart;

            out << APD_count *BCL << " "
                << APD90 << " "
                << APD75 << " "
                << APD50 << " "
                << APD30 << " "
                << Vmax  << " "
                << Vmin  << " "
                << dVdtmax << std::endl;

            if (APD_out_swhich == 0) {
                APD_out_end[0] = t - timeAPDstart;
                APD_out_swhich = 1;
            }
            else if (APD_out_swhich == 1) {
                APD_out_end[1] = t - timeAPDstart;
                APD_out_swhich = 0;
            }
        }
    }

    double dVdt = (Vm - Vm_prev) / dt;
    if (dVdt > dVdtmax) dVdtmax = dVdt;
    Istim_prev = Istim;
    Vm_prev = Vm;
}

/* measure the accumulations of calcium build up from the ICaL and RyR */
void APInfor::CalciumAccumulation(double t, double Istim, double BCL, double dt, double ICaL_con, double J_RyR) {
    if ((Istim >= 1e-3 || Istim <= -1e-3 ) && fabs(Istim_prev) < 1e-3) {

        out << t << " "
            << ICaL_in << " "
            << RyR_in  << std::endl;
        timeAPDstart = t;
        ICaL_in = 0.0;
        RyR_in = 0.0;
    }

    ICaL_in += ICaL_con * dt;
    RyR_in += J_RyR * dt;
    Istim_prev = Istim;
}

void APInfor::ReportAPD() {
    std::cout <<  "APD_out_end[0] = " <<  APD_out_end[0] << std::endl;
    std::cout <<  "APD_out_end[1] = " <<  APD_out_end[1] << std::endl;
}


void APInfor::ReportLastTwo() {
    std::cout /*<< BCL << "\t"*/
            << std::setprecision(9) << APD_out_end[0] << "\t"
            << APD75_out[0] << "\t"
            << APD50_out[0] << "\t"
            << APD30_out[0] << "\t"
            << Vmax_out[0] << "\t"
            << Vmin_out[0] << "\t"
            << dVdtmax_out[0] << "\t"
            << plateau_potential_out[0] << "\t"
            << INa_max_out[0] << "\t"
            << dVdt_Time_RP[0] << "\t"
            /* std::cout*/ /*<< BCL << "\t"*/
            << APD_out_end[1] << "\t"
            << APD75_out[1] << "\t"
            << APD50_out[1] << "\t"
            << APD30_out[1] << "\t"
            << Vmax_out[1] << "\t"
            << Vmin_out[1] << "\t"
            << dVdtmax_out[1] << "\t"
            << plateau_potential_out[1] << "\t"
            << INa_max_out[1] << "\t"
            << dVdt_Time_RP[1] << "\t"
            << V_m70_Time_RP[1]
            << std::endl;
}

void APInfor::ReportLastTwo(double BCL) {
    std::cout << std::setprecision(9) << BCL << "\t"
              << APD_out_end[0] << "\t"
              << APD75_out[0] << "\t"
              << APD50_out[0] << "\t"
              << APD30_out[0] << "\t"
              << APD20_out[0] << "\t"
              << Vmax_out[0] << "\t"
              << Vmin_out[0] << "\t"
              << dVdtmax_out[0] << "\t"
              << plateau_potential_out[0] << "\t"
              << INa_max_out[0] << "\t"
              << dVdt_Time_RP[0] << "\t"
              << V_m70_Time_RP[0] << "\t"
              << Cai_max_out[0] << "\t"
              << Cai_min_out[0] << "\t"
              /* std::cout*/ << BCL << "\t"
              << APD_out_end[1] << "\t"
              << APD75_out[1] << "\t"
              << APD50_out[1] << "\t"
              << APD30_out[1] << "\t"
              << APD20_out[1] << "\t"
              << Vmax_out[1] << "\t"
              << Vmin_out[1] << "\t"
              << dVdtmax_out[1] << "\t"
              << plateau_potential_out[1] << "\t"
              << INa_max_out[1] << "\t"
              << dVdt_Time_RP[1] << "\t"
              << V_m70_Time_RP[1] << "\t"
              << Cai_max_out[1] << "\t"
              << Cai_min_out[1]
              << std::endl;
}



void APInfor::ReportLast() {
    std::cout << APD90 << " "
              << APD75 << " "
              << APD50 << " "
              << APD30 << " "
              << Vmax  << " "
              << Vmin  << " "
              << dVdtmax << " "
              << plateau_potential << std::endl;
}

#endif