local_dir = 'C:\Users\timot\Desktop\trust_data';
cd(local_dir);

destination_path = 'C:\Users\timot\Desktop\trust_test';

id = dir;
id = id(~startsWith({id.name}, 'ex_'));

% %one sub
%id = id(strcmp({id.name}, '221603'));
%mm = id.name;
%trustbehavior(mm, destination_path, local_dir)

for i = 3:length(id)
    mm = id(i).name;
    trustbehavior(mm, destination_path, local_dir)
end
