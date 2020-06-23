clear
prot_rate=1;

for M = 0:11
for i = 0:7
		disp('rate = ');
		disp(prot_rate);
        disp([M, i]);
    for j = 1:600
        
        
        try
            measure_APD_From_mat_file(prot_rate, M,j,i)
        catch ME
        end
    end
end

end



%% functionname: function description
function [outputs] = measure_APD_From_mat_file(prot_rate, AF_model_ID, Cell_ID, Drug_ID)



file_AP = sprintf('AF_model_ID_%d/AP_Cell_%d_drug_%d.mat', AF_model_ID, Cell_ID, Drug_ID);
% file_AP = sprintf('AP_Cell_%d_drug_%d.mat',  Cell_ID, Drug_ID);

load (file_AP)

% Vm = y(:,39); Cion = y(:,38); Oion = y(:,66)+y(:,70);
% prot_rate=3;
Vm = Vmnew;
Cconc = Cconcnew;
Naconc =Naconcnew;

t = tnew;


% interplate the traces to a higher resolution data;
t_new_interval=min(tnew):0.05:max(tnew);

vq2 = interp1(tnew,Vmnew,t_new_interval,'spline');

t = t_new_interval';
Vm = vq2';

vq2 = interp1(tnew,Cconc,t_new_interval,'spline');
Cconc=vq2';

vq2 = interp1(tnew,Naconc,t_new_interval,'spline');
Naconc=vq2';




% Iion = currents(:,1); %dVm = conv2(Vm, [1;-1],'same');%currents(:,5);
period=1000/prot_rate;
dVm = (Vm - [Vm(1); Vm(1:end-1)])./(t - [t(1); t(1:end-1)]);%Vm - [Vm(1), Vm(1:end-1)];
% protocol_analysis_AP

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Analysis of Ca2+ properties (Vmax and APD) %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cconc = y(:,38); Naconc = y(:,34);
period = 1000/prot_rate;

% protocol_analysis_Ca

% AP_info.CA_t50 = t(t_in_index+index_max+CA_roi(1))-t(t_in_index+index);


AP_info1 = measure_AP_info(t, Vm,dVm,Cconc,period, Naconc, 1);
AP_info2 = measure_AP_info(t, Vm,dVm,Cconc,period, Naconc, 2);
% AF_model_ID=AF_model_ID;
CPU_ID=1;
drug_cond=Drug_ID;
model_ID=Cell_ID;
APfile = sprintf('Measure_from_file_AF_model_ID_%d_CPU_ID_%d_test_drug_%d.dat', AF_model_ID, CPU_ID, drug_cond);


OutputAPInfo(APfile, AP_info1, AP_info2, 'a+',model_ID);

t_in=t(end)-2*period-50; t_in_roi=find(t>t_in); t_in_index=t_in_roi(1)-1;


% file_AP = sprintf('AF_model_ID_%d_AP_Cell_%d_drug_%d.mat', AF_model_ID, model_ID, drug_cond);

% tnew = t(t_in_index:end);
% Vmnew = Vm(t_in_index:end);
% Cconcnew = Cconc(t_in_index:end);
% Naconcnew = Naconc(t_in_index:end);

% currents = calcCurrents(t(t_in_index:end),y(t_in_index:end, :),p); 
% figure,
% plot(tnew, currents(:,11))

% save(file_AP, 'tnew', 'Vmnew', 'Cconcnew', 'Naconcnew','currents');
end



% protocol_analysis_AP: last AP


function [AP_info] = measure_AP_info(t, Vm,dVm,Cconc, period, Na_conc, AP1)



t_in=t(end)-period; t_in_roi=find(t>t_in); t_in_index=t_in_roi(1)-1;

t_fin=t_in+period-5; t_fin_roi=find(t>t_fin); t_fin_index=t_fin_roi(1);
if AP1 == 1
    ;
else 
    t_in=t(end)-2*period; t_in_roi=find(t>t_in); t_in_index=t_in_roi(1)-1;

t_fin=t_in+period-5; t_fin_roi=find(t>t_fin); t_fin_index=t_fin_roi(1);
end







[val, index]=max(dVm(t_in_index:t_fin_index)); index=index-1; % max slope
dVm_max=val;
AP_info.dVdtmax = dVm_max;
[Vm_max, index_max]=max(Vm(t_in_index:t_fin_index)); index_max=index_max-1; % peak
Vm_min=min(Vm(t_in_index:t_fin_index)); % resting
AP_amp=Vm_max-Vm_min; 
% time_max = t_in(index_max);

APfraction_90=0.9;
APD_roi=find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction_90*AP_amp)-1;
AP_info.APD90 =t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index);
AP_info.APD90_max_amp =t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index_max);





% t_in=t(end)-period; t_in_roi=find(t>t_in); t_in_index=t_in_roi(1)-1;
% t_fin=t_in+period-5; t_fin_roi=find(t>t_fin); ...
%     t_fin_index=t_fin_roi(1);

plotfigure=0;
if(plotfigure == 1)

figure, set(gcf,'color','w')
subplot(2,1,1)
hold on,plot(t,dVm,'k',t(t_in_index:t_fin_index),dVm(t_in_index:t_fin_index),'r')
plot(t(t_in_index+index),dVm_max,'*')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('dEm (mV/ms)')
text(t(end)-period*2/3,dVm_max,['dEm max = ',num2str(dVm_max),' mV/ms']);
subplot(2,1,2)
hold on,plot(t,Vm,'k',t(t_in_index:t_fin_index),Vm(t_in_index:t_fin_index),'r')
plot(t(t_in_index+index_max),Vm_max,'b*')
plot([t(t_fin_index)-period t(t_fin_index)],[Vm_min Vm_min],'b:',[t(t_fin_index)-period t(t_fin_index)],[Vm_max-APfraction_90*AP_amp Vm_max-APfraction_90*AP_amp],'g:')
plot(t(t_in_index+index_max+APD_roi(1)),Vm(t_in_index+index_max+APD_roi(1)),'g*')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('Em (mV)'),xlabel('Time (ms)')
text(t(end)-period*2/3,Vm_max-AP_amp/3,['APD',num2str(100*APfraction_90),' = ',num2str(AP_info.APD90),' ms']);
% text(t(end),Vm_max-AP_amp/3, ['APD_a_m_p =', num2str(t(t_in_index+index_max+APD_roi(1)) -  t(t_in_index+index)),' ms']);


end


APDfraction_20=0.2;
APD_roi=find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APDfraction_20*AP_amp)-1;
AP_info.APD20=t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index);


APDfraction_75=0.75;
APD_roi=find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APDfraction_75*AP_amp)-1;
AP_info.APD75=t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index);

APDfraction_50=0.5;
APD_roi=find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APDfraction_50*AP_amp)-1;
AP_info.APD50=t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index);


AP_info.PLT20 = Vm_max-APDfraction_20*AP_amp;

AP_info.Vmax = Vm_max;
AP_info.Vmin = Vm_min;

AP_info.dVdtmax_time = t(t_in_index+index);

t_10_since_dvdt = find(t>t(t_in_index+index)+10);
t_10_since_dvdt_index=t_10_since_dvdt(1)-1;
t_50_since_dvdt = find(t>t(t_in_index+index)+50); 

t_50_since_dvdt_index=t_50_since_dvdt(1)-1;

AP_info.PLT = mean(Vm(t_10_since_dvdt_index:t_50_since_dvdt_index));




% Maximum calcium, minimum calcium, and amplitude
[CA_max, index_max] = max(Cconc(t_in_index:t_fin_index)); ...
    index_max=index_max-1; %peak
CA_min=min(Cconc(t_in_index:t_fin_index)); %resting
CA_amp = CA_max - CA_min;



CA_fraction = 0.5;
CA_roi=find(Cconc(t_in_index+index_max:t_fin_index)...
    <CA_max-CA_fraction*CA_amp)-1;
AP_info.CA_max = CA_max;
AP_info.CA_min = CA_min;

AP_info.Na_max = max( Na_conc(t_in_index:t_fin_index)); 
AP_info.Na_min = min( Na_conc(t_in_index:t_fin_index)); 


AP_info.CA_t50 = t(t_in_index+index_max+CA_roi(1))-t(t_in_index+index);






% APfraction=0.7
% APD_roi=find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
% APD=t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index)
% 
% APDfraction_50=0.5;
% APD_roi=find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APDfraction_50*AP_amp)-1;
% APD50 = t(t_in_index+index_max+APD_roi(1))-t(t_in_index+index);
% 
% figure, set(gcf,'color','w')
% subplot(2,1,1)
% hold on,plot(t,dVm,'k',t(t_in_index:t_fin_index),dVm(t_in_index:t_fin_index),'r')
% plot(t(t_in_index+index),dVm_max,'*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('dEm (mV/ms)')
% text(t(end)-period*2/3,dVm_max,['dEm max = ',num2str(dVm_max),' mV/ms']);
% subplot(2,1,2)
% hold on,plot(t,Vm,'k',t(t_in_index:t_fin_index),Vm(t_in_index:t_fin_index),'r')
% plot(t(t_in_index+index_max),Vm_max,'b*')
% plot([t(end)-period t(end)],[Vm_min Vm_min],'b:',[t(end)-period t(end)],[Vm_max-APDfraction_50*AP_amp Vm_max-APDfraction_50*AP_amp],'g:')
% plot(t(t_in_index+index_max+APD_roi(1)),Vm(t_in_index+index_max+APD_roi(1)),'g*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('Em (mV)'),xlabel('Time (ms)')
% text(t(end)-period*2/3,Vm_max-AP_amp/3,['APD',num2str(100*APDfraction_50),' = ',num2str(APD50),' ms']);
% 

% 
% figure, set(gcf,'color','w')
% subplot(2,1,1)
% hold on,plot(t,dVm,'k',t(t_in_index:t_fin_index),dVm(t_in_index:t_fin_index),'r')
% plot(t(t_in_index+index),dVm_max,'*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('dEm (mV/ms)')
% text(t(end)-period*2/3,dVm_max,['dEm max = ',num2str(dVm_max),' mV/ms']);
% subplot(2,1,2)
% hold on,plot(t,Vm,'k',t(t_in_index:t_fin_index),Vm(t_in_index:t_fin_index),'r')
% plot(t(t_in_index+index_max),Vm_max,'b*')
% plot([t(end)-period t(end)],[Vm_min Vm_min],'b:',[t(end)-period t(end)],[Vm_max-APDfraction_20*AP_amp Vm_max-APDfraction_20*AP_amp],'g:')
% plot(t(t_in_index+index_max+APD_roi(1)),Vm(t_in_index+index_max+APD_roi(1)),'g*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('Em (mV)'),xlabel('Time (ms)')
% text(t(end)-period*2/3,Vm_max-AP_amp/3,['APD',num2str(100*APDfraction_20),' = ',num2str(APD20),' ms']);

end




function OutputAPInfo(filename, APinfo, APinfo2, mode,model_ID)

file = fopen(filename, mode);

temp = '';
if nargin == 4
    temp = '';
elseif nargin == 5
    temp = model_ID;
end

fprintf(file, '%i ', temp);


fprintf(file, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f ', ... 
       APinfo.dVdtmax, ...
         APinfo.APD90, ...
 APinfo.APD90_max_amp, ...
         APinfo.APD20, ...
         APinfo.APD75, ...
         APinfo.APD50, ...
         APinfo.PLT20, ...
          APinfo.Vmax, ...
          APinfo.Vmin, ...
  APinfo.dVdtmax_time, ...
           APinfo.PLT, ...
        APinfo.CA_max, ...
        APinfo.CA_min, ...
        APinfo.CA_t50);


fprintf(file, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f ', ... 
       APinfo2.dVdtmax, ...
         APinfo2.APD90, ...
 APinfo2.APD90_max_amp, ...
         APinfo2.APD20, ...
         APinfo2.APD75, ...
         APinfo2.APD50, ...
         APinfo2.PLT20, ...
          APinfo2.Vmax, ...
          APinfo2.Vmin, ...
  APinfo2.dVdtmax_time, ...
           APinfo2.PLT, ...
        APinfo2.CA_max, ...
        APinfo2.CA_min, ...
        APinfo2.CA_t50);



% outputing Na info
fprintf(file, '%f %f %f %f\n', ... 
        APinfo.Na_max, ...
        APinfo.Na_min, ...
        APinfo2.Na_max, ...
        APinfo2.Na_min);



fclose(file);


end











