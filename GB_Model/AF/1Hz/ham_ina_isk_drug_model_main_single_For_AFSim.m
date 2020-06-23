function ham_ina_isk_drug_model_main_single_For_AFSim(CPU_ID, rate, AF, para,model_ID, AF_model_ID, AF_Remodl, AF_block, drug_cond, simulation_purpose)
%% Load initial conditions
% y0 contains the initial values of all state variables

%%%%%%%%%%%%%%%%%%%%%%%
% OLD MODELS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% w/ integrated IKur Markov model (w/o Isk)
% load yf_ha_orig_isk_1Hz % hh model
%load yf_ha_ikur_isk_1Hz % markov model

% load yf_ha_ikur_df_1Hz_AF_hh
% load yf_ha_ikur_df_3Hz_AF_hh

%%%%%%%%%%%%%%%%%%%
% TO RUN NO DRUG %%
%%%%%%%%%%%%%%%%%%%

% w/ integrated IKur Markov model (no Isk)

% load yf_ha_ikur_df_moreg_kall2_65_1Hz % markov model, 1 Hz
% load yf_ha_ikur_df_moreg_kall2_65_3Hz % markov model, 3 Hz



% default parameters here // 15:23:24, Wed, 02-May-2018, By Haibo
file_AP = sprintf('y_ini.mat'); % to get initial condition, run the model for 300 ms
duration = 300e3; % (ms)  run 300 s to get initial conditions


if  (~strcmp(simulation_purpose,'initial_condition')) % else, read initial conditons

	file_AP = sprintf('y_ini_ID_%d.mat', CPU_ID);
	duration = 500e3; % (ms)  run 500 s before measuring APD after drug.

end 

load (file_AP)  % get ini_conditions, either from y_ini.mat or initial condition run.






y0 = yfinal(1:86+3); % the new reduced model removes INa drug bound states and has less num of states; // 10:52:26, Tue, 19-November-2019, By Haibo


%% Simulation parameters

% Stimulation protocol parameters
prot_index = 1; % Selection of the protocol (-)
% 1) 'pace_cc'; 2) 'pace_cc_erp'; 3) 'tri_ap_clamp'; 4) 'ap_clamp';
% 5) 'v_hold'; 6) 'ssa_prot'; 7) 'ssi_prot'; 8) 'rec_prot';
% 9) 'v_step_tb'; 10) 'v_step'; 11) 'v_step_interval'; 12) 'rec_im_prot';
% 13) 'pace_cc_ead'; 14) 'pace_cc_dad'; 15) 'pace_cc_ead_sa';
% 16) 'pace_cc_ead_2beat'; 17) 'ap_clamp_w_step'; 18) 'ap_clamp_ramp'
prot_rate = rate; % PACING RATE (Hz) used only with 'pace_cc' and 'v_step'
prot_interval = 1000/3; % (ms) used only for "recovery" protocols (and ERP) (and 17)
prot_vm = yfinal(39); % (mV)

% Drug parameters
drug_index = 0; drug_conc = 0;               % Drug Free
%drug_index = 1; drug_conc = 10 * (1E-6);    % Ranolazine (M)
%drug_index = 2; drug_conc = 100 * (1E-6);    % Lidocaine (M)
%drug_index = 3; drug_conc = 10 * (1E-6);     % Flecainide (M)
%drug_index = 4; drug_conc = 0.3 * (1E-6);    % GS-967 (M)

% Other experimental conditions
exp_Temp = 310; % [K]
exp_Nao = 140; % [Na]o 130 mM in optimization, 140 otherwise
exp_ISO = 0; %; % (boolean)
exp_Ach = 0; %; % (boolean)   % no IKach

% IKur properties
% 0 drug free; 1 drug-CLOSED; 2 drug-OPEN; 3 drug-INACTIVATED; 4 drug-OPEN
% & INACTIVE
drug_index_ikur = 0; %(-)
Constant_Block = 0; % 1 if consant 50% block, 0 if not
Complete_Block = 0;
% drug_conc_kur = drug_conc_C_nr(1); % mM
drug_conc_kur = 0; % mM
kon = 0.00001; % 1/ms
koff = 0.00001;  % 1/ms

% Parameter array for passing nondefault conditions
prot_par = [prot_index prot_rate prot_interval prot_vm];    % 1 2 3 4
drug_par = [drug_index drug_conc];                          % 5 6
exp_par = [exp_Temp exp_Nao exp_ISO exp_Ach];               % 7 8 9 10
ikur_par = [drug_index_ikur Constant_Block drug_conc_kur kon koff Complete_Block]; % 11 12 13 14 15 16
% AF = 0;


 % 1 = AF, 0 = no AF
p.option = [prot_par, drug_par, exp_par, ikur_par, AF]; %1


p.ptimeVC=0; 
p.emVC_1=0;  
p.emVC_2=0; 
p.EAD_Vm=0; 
p.EAD_t=0; 
p.altTime=0; 
p.altCaj=0; 
p.altCasls=0;
% p.ScaleFactors = ones(14, 1);


S_GNa = 1;
S_GNaB = 1;
S_GNaK = 1;
S_GClCa = 1;
S_GCaL = 1;
S_GCaB = 1;
S_GNCX = 1;
S_GCaP = 1;
S_GKs = 1;
S_GKAch = 1;
S_Gto = 1;
S_GKur = 1;
S_GK1 = 1;
S_GKr = 1;
S_GClB = 1;
S_GSK = 1;
S_GK2P = 1;
% p.ScaleFactors = [S_GNa,S_GClCa, S_GCaL, S_Gto, S_GKur,S_GKr,S_GKs,  S_GK1, S_GNCX,S_GNaK, S_GSK, S_GK2P , S_GNaB, S_GCaB, S_GCaP,   S_GClB ];  % , S_GKAch,  not updated 
% p.data = [ptimeVC emVC_1 emVC_2 EAD_Vm EAD_t altTime altCaj altCasls];  % 

p.ScaleFactors = para;


p.S_GKur_AF = AF_Remodl.S_GKur_AF * (1 - AF_block.S_GKur_AF);%0.5;
p.S_GSK_AF  = AF_Remodl.S_GSK_AF* (1 - AF_block.S_GSK_AF);%2.0;
p.S_GK2P_AF = AF_Remodl.S_GK2P_AF* (1 - AF_block.S_GK2P_AF);%3.0;


% global timeVC emVC_1 emVC_2 EAD_Vm EAD_t
% global altTime altCaj altCasl

% Simulation duration
% duration = 2e3/prot_rate; % (ms)


%% Single Run Simulation
% y0 = [y0,0,0,0];

tic
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','off'); 
[t,y] = ode15s(@ham_ina_isk_drug_model_reduce,tspan,y0,options,p);
yfinal = y(end,:);
toc




%% well, only compute these values when doing things other than initial conditions purpose.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of AP properties (Vmax and APD) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  strcmp(simulation_purpose,'initial_condition') % save yfinal after calculating initial_conditions.

	file_AP = sprintf('y_ini_ID_%d.mat', CPU_ID);
	save(file_AP, 'yfinal');

else % when not calculating initial_conditions.
	if (~isreal(y))
		y = real(y);
	end
	Vm = y(:,39); Cion = y(:,38); Oion = y(:,66)+y(:,70);
	% Iion = currents(:,1); %dVm = conv2(Vm, [1;-1],'same');%currents(:,5);
	period=1000/prot_rate;
	dVm = (Vm - [Vm(1); Vm(1:end-1)])./(t - [t(1); t(1:end-1)]);%Vm - [Vm(1), Vm(1:end-1)];
	% protocol_analysis_AP

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% % Analysis of Ca2+ properties (Vmax and APD) %%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Cconc = y(:,38); Naconc = y(:,34);
	period = 1000/prot_rate;

	% protocol_analysis_Ca

	% AP_info.CA_t50 = t(t_in_index+index_max+CA_roi(1))-t(t_in_index+index);


	AP_info1 = measure_AP_info(t, Vm,dVm,Cconc,period, Naconc, 1);
	AP_info2 = measure_AP_info(t, Vm,dVm,Cconc,period, Naconc, 2);


	APfile = sprintf('AF_model_ID_%d/CPU_ID_%d_test_drug_%d.dat', AF_model_ID, CPU_ID, drug_cond);


	OutputAPInfo(APfile, AP_info1, AP_info2, 'a+',model_ID);

	t_in=t(end)-2*period-50; t_in_roi=find(t>t_in); t_in_index=t_in_roi(1)-1;


	file_AP = sprintf('AF_model_ID_%d/AP_Cell_%d_drug_%d.mat', AF_model_ID, model_ID, drug_cond);

	tnew = t(t_in_index:end);
	Vmnew = Vm(t_in_index:end);
	Cconcnew = Cconc(t_in_index:end);
	Naconcnew = Naconc(t_in_index:end);

	currents = calcCurrents(t(t_in_index:end),y(t_in_index:end, :),p); 
	% figure,
	% plot(tnew, currents(:,11))

	save(file_AP, 'tnew', 'Vmnew', 'Cconcnew', 'Naconcnew','currents');


end 
toc

end