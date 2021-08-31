function [vo, poolobj, behavfiles] = trust_setup_environment(vo)
%this function sets up the paths for sceptic, VBA, and expected behavior files based on the user's environment
%it also handles setup of the parpool depending on the user

%os = computer;
[~, me] = system('whoami');
me = strtrim(me);

[~, host]=system('hostname');
host = strtrim(host);

is_aci = contains(host, 'aci.ics.psu.edu');
is_longleaf = contains(host, 'unc.edu');

is_alex = strcmp(me,'Alex')==1 || strcmp(me,'dombax')==1;
is_jiazhouchen=strcmp(me,'jiazhouchen')==1;
%is_jiazhou_t=strcmp(me,'jiazhou')==1;
is_mnh=(contains(me, 'mnh5174') || contains(me, 'michael')) || contains(me, 'mnhallq');
is_alison_mbp = strcmp(me, 'alisonmarie526') == 1;
is_alison_mac = strcmp(me, 'ams939') == 1;
is_daniel = strcmp(me,'dpp5430') == 1;
is_tim = strcmp(me, 'timallen') == 1 || strcmp(me, 'timot') == 1 || strcmp(me, 'timothya') == 1;

%% set environment and define file locations

if is_tim
    if is_longleaf
        vo.ta_vba_path='/proj/mnhallqlab/users/tim/mh_vba'; %MH helper scripts
        vo.vba_path = '/proj/mnhallqlab/lab_resources/VBA-toolbox';
        vo.trust_repo = '/proj/mnhallqlab/users/tim/GitHub/Social_Trust'; 
    elseif is_mac
        vo.ta_vba_path='/proj/mnhallqlab/users/tim/mh_vba'; %MH helper scripts
        vo.vba_path = '/proj/mnhallqlab/lab_resources/VBA-toolbox';
        vo.trust_repo = '/proj/mnhallqlab/users/michael/Social_Trust'; 
    else
        vo.ta_vba_path='/proj/mnhallqlab/users/tim/mh_vba'; %MH helper scripts
        vo.vba_path = '/proj/mnhallqlab/lab_resources/VBA-toolbox';
        vo.trust_repo = '/proj/mnhallqlab/users/michael/Social_Trust'; 
    end
elseif is_alex
    vo.vba_path='~/code/VBA-toolbox';
    vo.trust_repo='~/code/Social_Trust';
    vo.mh_vba_path='~/code/mh_vba'; %MH helper scripts
    vo.output_dir = '~/Box/skinner/projects_analyses/Trust/mfx_analyses';
    vo.data_directory = [vo.trust_repo,'/data/',vo.dataset,'/mat'];
elseif is_jiazhouchen
    vo.vba_path='/Volumes/bek/vbatoolbox';
    vo.boxdir = '~/Box';
    vo.output_dir = '~/Box/skinner/projects_analyses/Trust/mfx_analyses';
elseif is_mnh
    if is_aci
        vo.mh_vba_path='/gpfs/group/mnh5174/default/Michael/mh_vba'; %MH helper scripts
        vo.vba_path = '/storage/home/mnh5174/MATLAB/VBA-toolbox';
        vo.trust_repo = '/gpfs/group/mnh5174/default/Michael/Social_Trust';
    elseif is_longleaf
        vo.mh_vba_path='/proj/mnhallqlab/users/michael/mh_vba'; %MH helper scripts
        vo.vba_path = '/proj/mnhallqlab/lab_resources/VBA-toolbox';
        vo.trust_repo = '/proj/mnhallqlab/users/michael/Social_Trust';      
    else
        %local mac installation
        vo.mh_vba_path='/Users/michael/Data_Analysis/mh_vba'; %MH helper scripts
        vo.vba_path = '~/Documents/MATLAB/VBA-toolbox';
        vo.trust_repo='/Users/michael/Data_Analysis/Social_Trust';
    end
       
    vo.data_directory=[vo.trust_repo, '/data/', vo.dataset, '/mat'];
    vo.output_dir=[vo.trust_repo, '/output/vba'];
    %vo.output_dir=['/Users/michael/trust_tmp/output/vba'];
elseif is_alison_mbp
    vo.vba_path='/Users/alisonmarie526/Downloads/VBA-toolbox-master-2';
    vo.trust_repo='/Users/alisonmarie526/github_dirs/Social_Trust';
    vo.output_dir = [vo.trust_repo, '/output/vba'];
    vo.data_directory = [vo.trust_repo,'/data/',vo.dataset,'/mat'];
    mmrescan_dir = '/Users/alisonmarie526/github_dirs/Social_Trust';
elseif is_alison_mac
    vo.vba_path='/Users/ams939/Downloads/VBA-toolbox-master-2';
    vo.trust_repo = '~/Box Sync/DEPENd/Projects/Social_Trust';
    vo.output_dir = [vo.trust_repo, '/output/vba'];
    mmrescan_dir = '/Users/ams939/Box Sync/DEPENd/Projects/Social_Trust';
elseif is_daniel
    if is_aci
        vo.mh_vba_path='/gpfs/group/mnh5174/default/Michael/mh_vba'; %MH helper scripts
        vo.vba_path='/storage/home/dpp5430/VBA-toolbox';
        vo.trust_repo = '/storage/home/dpp5430/Social_Trust';
        vo.output_dir = '/storage/home/dpp5430/ST_output/special_test';
    elseif is_longleaf
        vo.mh_vba_path='/proj/mnhallqlab/users/michael/mh_vba'; %MH helper scripts
        vo.vba_path = '/storage/home/mnh5174/MATLAB/VBA-toolbox';
        vo.trust_repo = '/proj/mnhallqlab/users/michael/Social_Trust';      
    else
        %vo.mh_vba_path=PUT ON MAC
        vo.vba_path='/Users/dpp5430/Desktop/VBA-toolbox';
        vo.trust_repo = '/Users/dpp5430/Desktop/workingST';
        vo.output_dir = '/Users/dpp5430/Desktop/final_model_contestants';
    end 
    vo.data_directory = [vo.trust_repo,'/data/',vo.dataset,'/mat'];
else
    vo.output_dir = [vo.trust_repo, '/output/vba_out/', vo.dataset, '/mfx/', vo.model];
    if ~exist(vo.output_dir, 'dir'), mkdir(vo.output_dir); end
end

addpath(vo.mh_vba_path); %add MH VBA helpers (e.g., extract_group_statistics)
addpath(genpath_safe(vo.vba_path)); %add VBA functions

if ~isfield(vo, 'dataset'), vo.dataset='mmrescan'; end %default
if ~isfield(vo, 'subj_subset'), vo.subj_subset=''; end %fit all by default


if strcmpi(vo.dataset, 'mmrescan') 
    % identify subject behavior files
    demographics_file = readtable([vo.trust_repo, '/Rescanning_Project_Scan_Log.csv']);
    %demographics_file = removevars(demographics_file, {'Notes', 'Confidence'});
    matching_file = readtable([vo.trust_repo, '/subsetting_pts/mmrescan_matchingresults_rist_jan2019.csv']);
    %matching_file = removevars(matching_file, {'Var1', 'matchid'});
    if strcmpi(vo.subj_subset, 'normative')
        %dataset containing all controls
        controls=strcmpi(demographics_file.Dx_Group, 'HC');
        behavfiles=strcat(vo.data_directory, '/trust', demographics_file.MR_Center_ID(controls), '.mat');
    elseif strcmpi(vo.subj_subset, 'normative62')
        %all controls except two extreme cases (who shared almost always)
        controls=strcmpi(demographics_file.Dx_Group, 'HC');
        ids=demographics_file.MR_Center_ID(controls);
        filter_vec = ~ismember(ids, {'111jv_14Dec2016', '11252_20161118'});
        ids = ids(filter_vec);
        behavfiles=strcat(vo.data_directory, '/trust', ids, '.mat');        
    elseif strcmpi(vo.subj_subset, 'matched')
        behavfiles=strcat(vo.data_directory, '/trust', matching_file.MR_Center_ID, '.mat');
    elseif strcmpi(vo.subj_subset, 'specc')
        %dataset containing both BPD and matched control participants
        specc=contains(demographics_file.Study, 'SPECC');
        behavfiles=strcat(vo.data_directory, '/trust', demographics_file.MR_Center_ID(specc), '.mat');
    else
        fprintf('Fitting all subjects because subj_subset not used/matched\n');
        behavfiles=strcat(vo.data_directory, '/trust', demographics_file.MR_Center_ID, '.mat');
    end
    
    %verify that expected files exist, and remove nonexistent files with warning
    expect_exist=cellfun(@(x) exist(x, 'file'), behavfiles);
    if ~all(expect_exist)
        miss_files=(expect_exist == 0);
        fprintf('\nThe following expected files are missing:\n');
        fprintf('   %s\n', behavfiles{miss_files});
        fprintf('\nRemoving these from fitting for the time being, but make sure we resolve discrepancies!\n\n');
        behavfiles=behavfiles(~miss_files);        
    end   
    
elseif strcmpi(vo.dataset, 'bsocial')
    vo.data_directory=fullfile(vo.boxdir,'skinner/data/eprime/bpd_trust');
    addpath import;
    IDdirs=dir(vo.data_directory);
    IDs={IDdirs.name};
    IDs=IDs(~cellfun(@(x) strcmpi(x,'.') | strcmpi(x,'..'),IDs));
    allfilepath = []; allb=[];
    for id = IDs
        IDtext=sprintf('ID%s', id{:});
        [allb.(IDtext),allfilepath.(IDtext)] = trustbehavior(vo.data_directory,id{:},0);
    end
    vo.allb=allb;
    vo.n_trials=144;
    behavfiles=struct2cell(allfilepath);
    behavfiles=behavfiles(~cellfun(@isempty,behavfiles));
elseif strcmpi(vo.dataset,'learn')
    vo.data_directory=fullfile(vo.boxdir,'skinner/data/eprime/scan_trust');
    addpath import;
    IDs=[{'35780'}, {'46069'}, {'202200'}, {'202278'}, {'203264'}, {'203803'}, {'204015'}, {'206270'}, {'207224'}, {'208510'}, {'210290'}, {'210374'}, {'210701'}, {'211038'}, {'211101'}, {'211253'}, {'212019'}, {'212020'}, {'212794'}, {'213163'}, {'213704'}, {'213868'}, {'214082'}, {'215088'}, {'215537'}, {'215622'}, {'216011'}, {'216032'}, {'217008'}, {'217048'}, {'217173'}, {'217293'}, {'217988'}, {'218219'}, {'218572'}, {'218714'}, {'219044'}, {'219089'}, {'219392'}, {'219619'}, {'219658'}, {'219676'}, {'219772'}, {'219809'}, {'219886'}, {'219944'}, {'219949'}, {'219954'}, {'220043'}, {'220051'}, {'220064'}, {'220104'}, {'220152'}, {'220184'}, {'220186'}, {'220244'}, {'220299'}, {'220399'}, {'220451'}, {'220477'}, {'220513'}, {'220523'}, {'220531'}, {'220553'}, {'220566'}, {'220597'}, {'220678'}, {'220684'}, {'220740'}, {'220758'}, {'220789'}, {'220796'}, {'220865'}, {'220889'}, {'220913'}, {'220926'}, {'220927'}, {'220947'}, {'220963'}, {'220976'}, {'220989'}, {'221018'}, {'221036'}, {'221050'}, {'221085'}, {'221099'}, {'221140'}, {'221164'}, {'221181'}, {'221183'}, {'221206'}, {'221212'}, {'221226'}, {'221261'}, {'221273'}, {'221292'}, {'221295'}, {'221298'}, {'881075'}, {'881105'}, {'881106'}, {'881132'}, {'881224'}, {'881230'}];
    allfilepath = []; allb=[];
    for id = IDs
        IDtext=sprintf('ID%s', id{:});
        [allb.(IDtext),allfilepath.(IDtext)] = trustbehavior(vo.data_directory,id{:},0);
    end
    vo.allb=allb;
    vo.n_trials=192;
    behavfiles=struct2cell(allfilepath);
    behavfiles=behavfiles(~cellfun(@isempty,behavfiles));
end


%% setup parallel parameters
if ~isfield(vo, 'do_parallel')
    vo.ncpus = 1;
    vo.do_parallel = false;
elseif ~vo.do_parallel
    vo.ncpus = 1;
end

%setup per-user defaults
if is_alex
    if strcmp(me,'dombax')==1
        vo.ncpus = 12;
        fprintf('defaulting to 12 cpus on Thorndike \n');
    else
        vo.ncpus=10;
        fprintf('defaulting to 10 cpus on iMac Pro \n');
    end
elseif is_jiazhouchen
    vo.ncpus=4;
elseif is_daniel
    if is_aci
        %vo.ncpus=8;
        vo.ncpus=16;
    end
end

ncpus=getenv('matlab_cpus');
if ~strcmpi(ncpus, '')
    vo.ncpus = str2double(ncpus);
    fprintf('Using %d cpus based on matlab_cpus environment variables\n', vo.ncpus);
end

if vo.ncpus > 1, vo.do_parallel=true; end

if vo.do_parallel    
    delete(gcp('nocreate')); %stop any pool that is running
    poolobj=parpool('local',vo.ncpus); %just use shared pool for now since it seems not to matter (no collisions)
else
    poolobj=[];
end

end
