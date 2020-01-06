function b = asterisk(subdir)
%subdir = 'f_trust_Qlearn_counter_hybrid';
%adding asterisk to existing .dat files
%subdir=subdir{:};
%data_dir_str= sprintf('/Users/polinavanyukov/Regs/clinical_policy_values/01-Aug-2019/');
%data_dump_str = data_dir_str;

%data_dir_str=['/Users/polinavanyukov/Regs/clinical_policy_PEsValues_FixedPs/',date];
data_dir_str=['/Users/polinavanyukov/Regs/clinical_policy_PEsValues/',date];
data_dump_str =['/Users/polinavanyukov/Regs/asterisked_regs/',date,'/'];


if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('*.dat');
num_of_subjects = length(files);


for index = 1:num_of_subjects
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    x = load(filename);
    block1=num2cell(x(1:48,:));
    block2=num2cell(x(49:96,:));
    block3=num2cell(x(97:144,:));
    block4=num2cell(x(145:192,:));
    ast = {'*', '*', '*'};
    c = [block1; ast; block2; ast; block3; ast; block4];
    %writetable(cell2table(c), [data_dump_str filename],'Delimiter','\t');
    dlmcell([data_dump_str filename],c,'\t');
end


%Copy each models PE's to a date dir somewhere to consolidate it for easy
%movement to Thorndike
%to_thorndike_dir='E:\data\trust\regs\last_pes\';
%copyfile([data_dump_str '*.dat'],[to_thorndike_dir date filesep subdir filesep]);

return

