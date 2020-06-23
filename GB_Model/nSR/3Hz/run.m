
para = load('para_log.dat');



parfor ind = 1:600

ind

p = para(ind, :);

t=getCurrentTask();

ID=t.ID;
try
ham_ina_isk_drug_model_main_single(ID, 3, 0, p, ind);

catch ME
end
end