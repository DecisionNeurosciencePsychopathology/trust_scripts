function b = align_PEs_decision()

data_dir_PEs= '/Users/polinavanyukov/Regs/hc_policy_PEs_nonasterisk';
data_dir_timings = '/Users/polinavanyukov/Regs/hc_decision_nonasterisk';
data_dump_str = '/Users/polinavanyukov/Regs/PEs_aligned_decision';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_PEs)
files1 = dir('trust*.dat');
files2 = dir([data_dir_timings, '/*decision_Times.dat']);
num_of_subjects = length(files1);

for index=1:num_of_subjects
    filename1 = files1(index).name;
    filename2 = files2(index).name;
    id = filename1(isstrprop(filename1,'digit'));
    PEs = load(filename1);
    decision_times = load([data_dir_timings,'/',filename2]);
    newPEs = PEs;
    newPEs(:,1) = decision_times(:,1);
    data_dump_str_id = [data_dump_str,'/',num2str(id),'PEs_at_decision'];
    gdlmwrite(data_dump_str_id,newPEs,'\t');
end

end