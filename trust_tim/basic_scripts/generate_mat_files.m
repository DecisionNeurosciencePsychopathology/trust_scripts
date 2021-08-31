local_dir = 'C:\Users\timot\Box\skinner\personal_folders\Tim\trust_bsocial\data';
cd(local_dir);

mkdir 'C:\Users\timot\Box\skinner\personal_folders\Tim\trust_bsocial\data\processed'

destination_path = 'C:\Users\timot\Box\skinner\personal_folders\Tim\trust_bsocial\data\processed';
id = dir;
id = id(~startsWith({id.name}, 'ex_')); %exclude a couple folders that do not contain relevant subject data

% %one sub
%id = id(strcmp({id.name}, '221603'));
%mm = id.name;
%trustbehavior(mm, destination_path, local_dir)

for i = 3:length(id)
    mm = id(i).name;
    trustbehavior(mm, destination_path, local_dir)
end

cd(destination_path)
mkdir group_data

trust_group %process into a single group data file
