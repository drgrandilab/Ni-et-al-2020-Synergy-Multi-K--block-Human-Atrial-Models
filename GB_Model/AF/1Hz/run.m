
%% functionname: function description
function run(AF_model_ID)
    try
	pool=parpool(28);  % to check if give identical results when use single CPU
    catch
    ;    
    end
    
	para = load('para_log.dat');

	parfor i = 1:600
		p = para(i, :);

		t=getCurrentTask();

		ID=t.ID;
		fprintf('currently working on AF_Model %i, Cell model ID %i \n', AF_model_ID, i);
		% ham_ina_isk_drug_model_main_single(ID, 1, 0, p, i);



		% compute_single_para(1, 1, 1, ones(16,1), 1, 0)
		try

		compute_single_para(ID, 1, 1, p, i,AF_model_ID)
	catch ME
		disp(ME)
	end
	end

end