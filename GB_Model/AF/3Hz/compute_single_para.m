%% compute_single_para: function description
function [outputs] = compute_single_para(CPU_ID, rate, AF, para, model_ID,AF_Model_ID)


AF_remodelling =  dlmread('AF_models.dat');
Drug_cond =  dlmread('block_para.dat');



try
    
    newdir = sprintf('AF_model_ID_%d/', AF_Model_ID);
    if exist(newdir, 'dir') ~= 7
        mkdir(newdir);
    end
catch
    ;
end

% parameters in this order, identical to C++ code
% GKur
% GSK 
% GK2P


AF_Remodl.S_GKur_AF  = AF_remodelling(AF_Model_ID+1,1);  %;
AF_Remodl.S_GSK_AF = AF_remodelling(AF_Model_ID+1,2) ; %;
AF_Remodl.S_GK2P_AF = AF_remodelling(AF_Model_ID+1,3);  %;

AF_block.S_GKur_AF = Drug_cond(0+1, 1);%;
AF_block.S_GSK_AF = Drug_cond(0+1, 2);%;
AF_block.S_GK2P_AF = Drug_cond(0+1, 3);%;

simulation_purpose  = 'initial_condition';

ham_ina_isk_drug_model_main_single_For_AFSim(CPU_ID, rate, AF, para,model_ID, AF_Model_ID, AF_Remodl, AF_block, 0, simulation_purpose)


for i = 0:7


   fprintf('drug %i Cell ID %i AF Model ID %i\n', i, model_ID, AF_Model_ID);

AF_block.S_GKur_AF = Drug_cond(i+1, 1);
AF_block.S_GSK_AF = Drug_cond(i+1, 2);
AF_block.S_GK2P_AF = Drug_cond(i+1, 3);

simulation_purpose = 'drug';
ham_ina_isk_drug_model_main_single_For_AFSim(CPU_ID, rate, AF, para,model_ID, AF_Model_ID, AF_Remodl, AF_block, i, simulation_purpose)

end




end
