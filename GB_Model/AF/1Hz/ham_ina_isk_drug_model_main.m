%% ham_ina_isk_drug_model_main

clear all;
close all;
clc;

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
load y_ini.mat
% load yf_ha_ikur_df_1Hz_constantblock_299 % constant 50% pre block 1 Hz
% load yf_ha_ikur_df_3Hz_constantblock_300 % constant 50% pore IKur block, 3 Hz
% load yf_ha_ikur_df_1Hz_completeblock_300 % complete IKur block, 1 Hz
% load yf_ha_ikur_df_3Hz_completeblock_300 % complete IKur block, 3 Hz

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AF model - 3 Hz %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load yf_ha_ikur_df_3Hz_AF % AF, 3Hz, drug-free - steady-state took 900e3 ms
% load yf_ha_ikur_constant50_3Hz_AF %AF, constant50 block in IKur
% load yf_ha_ikur_completeblock_3Hz_AF %AF, complete block in IKur

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AF model - 1 Hz %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load yf_ha_ikur_df_1Hz_AF % AF, 1Hz, drug-free - also used 900e3 ms for steady-state
% load yf_ha_ikur_constant50_1Hz_AF %AF, constant50 block in IKur
% load yf_ha_ikur_completeblock_1Hz_AF %AF, complete block in IKur

y0 = yfinal;


%% Simulation parameters

% Stimulation protocol parameters
prot_index = 1; % Selection of the protocol (-)
% 1) 'pace_cc'; 2) 'pace_cc_erp'; 3) 'tri_ap_clamp'; 4) 'ap_clamp';
% 5) 'v_hold'; 6) 'ssa_prot'; 7) 'ssi_prot'; 8) 'rec_prot';
% 9) 'v_step_tb'; 10) 'v_step'; 11) 'v_step_interval'; 12) 'rec_im_prot';
% 13) 'pace_cc_ead'; 14) 'pace_cc_dad'; 15) 'pace_cc_ead_sa';
% 16) 'pace_cc_ead_2beat'; 17) 'ap_clamp_w_step'; 18) 'ap_clamp_ramp'
prot_rate = 1; % PACING RATE (Hz) used only with 'pace_cc' and 'v_step'
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
AF = 0;


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
S_GKur = 0.5;
S_GK1 = 1;
S_GKr = 1;
S_GClB = 1;
S_GSK = 1;
S_GK2P = 1;
p.ScaleFactors = [S_GNa,S_GClCa, S_GCaL, S_Gto, S_GKur,S_GKr,S_GKs,  S_GK1, S_GNCX,S_GNaK, S_GSK, S_GK2P , S_GNaB, S_GCaB, S_GCaP,   S_GClB ];  % , S_GKAch,  not updated 
% p.data = [ptimeVC emVC_1 emVC_2 EAD_Vm EAD_t altTime altCaj altCasls];  % 


% global timeVC emVC_1 emVC_2 EAD_Vm EAD_t
% global altTime altCaj altCasl

% Simulation duration
% duration = 2e3/prot_rate; % (ms)
duration = 100e3/prot_rate; % (ms)

%% Single Run Simulation
% y0 = [y0,0,0,0];

tic
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','on'); 
[t,y] = ode15s(@ham_ina_isk_drug_model,tspan,y0,options,p);
yfinal = y(end,:);
toc


% figure,plot(t/1000,y(:,34),'b')
%break
%% calcCurrents
% this function extracts the values of simulation outputs that are not state
% variables (e.g., currents) and we want to plot (selected in ham_ina_isk_drug_model)



%% Plot simulation outputs
% 
figure % membrane potential and ion concentrations are state variables
set(gcf,'color','w','Position',[100 100 600 800])
subplot(4,1,1); hold on, plot(t/1000,y(:,39),'k')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('E_m (mV)')
subplot(4,1,2); hold on, plot(t/1000,y(:,38)*1000,'b')
    plot(t/1000,y(:,36)*1000,'g',t/1000,y(:,37)*1000,'c')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Ca^2^+]_i (\muM)')
subplot(4,1,3); hold on, plot(t/1000,y(:,31)*1000,'b')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Ca^2^+]_S_R (\muM)')
subplot(4,1,4); hold on, plot(t/1000,y(:,34),'b')
    plot(t/1000,y(:,32),'g',t/1000,y(:,33),'c')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Na^+]_i (mM)')
xlabel('Time (s)')

% figure
% set(gcf, 'color', 'w')
% plot(t, y(:,105), 'r', t, 1 - y(:,101) - y(:,102) - y(:,103) - y(:,104) - y(:,105) - y(:,106) ...
%     - y(:,107) - y(:,108) - y(:,109) - y(:,110) - y(:,111), 'b', ...
%     t, y(:,101) + y(:,102) + y(:,103) + y(:,104), 'g', t, y(:,106) + y(:,107) + y(:,108) + y(:,109) ...
%     + y(:,110) + y(:,111), 'k');
% set(gca,'box','off','tickdir','out','fontsize',30)
% xlabel('Time (ms)')
% ylabel('State Occupancy')
% ylim([0 1])


% figure % currents are not state variables
% set(gcf,'color','w','Position',[100+700 100 600 800])
% subplot(4,1,1); hold on, plot(t/1000,y(:,39),'k')
% set(gca,'box','off','tickdir','out','fontsize',30)
% ylabel('E_m (mV)')
% subplot(4,1,2); hold on, plot(t/1000,y(:,38)*1000,'b')
% set(gca,'box','off','tickdir','out','fontsize',30)
% ylabel('[Ca^2^+]_i (\muM)')
% subplot(4,1,3); hold on, plot(t/1000,currents(:,3),'b',t/1000,currents(:,4),'g',t/1000,currents(:,1),'r')
% set(gca,'box','off','tickdir','out','fontsize',30)
% legend('I_C_a','I_N_C_X','I_N_a'), ylabel('Current (A/F)')
% subplot(4,1,4); hold on, plot(t/1000,currents(:,6),'k');
% set(gca,'box','off','tickdir','out','fontsize',30);
% ylabel('I_S_K (A/F)'); xlabel('Time (s)');

% figure % currents are not state variables
% set(gcf,'color','w')%,'Position',[100+700 100 600 800])
% subplot(2,1,1); hold on, plot(t,y(:,39),'k')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('E_m (mV)')
% subplot(2,1,2); hold on, plot(t,currents(:,8),'k')%, t, currents(:,9), 'b')
% % subplot(2,1,2); hold on, plot(t, currents(:,9), 'b')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('I_K_u_r (A/F)')
% xlabel('Time (ms)');

% figure(1) % currents are not state variables
% set(gcf,'color','w')%,'Position',[100+700 100 600 800])
% subplot(3,1,1); hold on, plot(t,y(:,39),'-')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('E_m (mV)')
% subplot(3,1,2); hold on, plot(t,currents(:,8),'-')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('I_K_u_r (A/F)')
% subplot(3,1,3); hold on, plot(t,y(:,38)*1000,'-')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('[Ca^2^+]_i (\muM)')
% xlabel('Time (ms)');

% figure(1) % currents are not state variables
% set(gcf,'color','w')%,'Position',[100+700 100 600 800])
% subplot(4,1,1); hold on, plot(t,y(:,39),'-')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('E_m (mV)')
% subplot(4,1,2); hold on, plot(t,currents(:,8),'-')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('I_K_u_r (A/F)')
% subplot(4,1,3); hold on, plot(t,y(:,38)*1000,'-')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('[Ca^2^+]_i (\muM)')
% xlabel('Time (ms)');
% subplot(4,1,4); hold on, plot(t, y(:,105), 'b', t, 1 - y(:,101) - y(:,102) - y(:,103) - y(:,104) - y(:,105), 'g', ...
%     t, y(:,101) + y(:,102) + y(:,103) + y(:,104), 'r');
% set(gca,'box','off','tickdir','out','fontsize',30)
% xlabel('Time (ms)')
% ylabel('State Occupancy')
% legend('Open', 'Inactive', 'Closed')  
% ylim([0 1])

% figure(1) % currents are not state variables
% set(gcf,'color','w')%,'Position',[100+700 100 600 800])
% subplot(3,1,1); hold on, plot(t,y(:,39),'-k')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('E_m (mV)')
% subplot(3,1,2); hold on, plot(t,currents(:,8),'-k')
% set(gca,'box','off','tickdir','out','fontsize',22)
% ylabel('I_K_u_r (A/F)')
% xlabel('Time (ms)');
% subplot(3,1,3); hold on, plot(t, y(:,105), 'b', t, 1 - y(:,101) - y(:,102) - y(:,103) - y(:,104) - y(:,105), 'g', ...
%     t, y(:,101) + y(:,102) + y(:,103) + y(:,104), 'r');
% set(gca,'box','off','tickdir','out','fontsize',22)
% xlabel('Time (ms)')
% ylabel('State Occupancy')
% legend('Open', 'Inactive', 'Closed')  
% ylim([0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of AP properties (Vmax and APD) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


AP_info1 = measure_AP_info(t, Vm,dVm,Cconc,period,1);
AP_info2 = measure_AP_info(t, Vm,dVm,Cconc,period,2);


OutputAPInfo('test', AP_info1, AP_info2, 'a+');

t_in=t(end)-2*period-50; t_in_roi=find(t>t_in); t_in_index=t_in_roi(1)-1;



file_AP = 'AP.mat';

tnew = t(t_in_index:end);
Vmnew = Vm(t_in_index:end);
Cconcnew = Cconc(t_in_index:end);

currents = calcCurrents(t(t_in_index:end),y(t_in_index:end, :),p); 
figure,
plot(tnew, currents(:,11))

save(file_AP, 'tnew', 'Vmnew', 'Cconcnew', 'currents');
toc




