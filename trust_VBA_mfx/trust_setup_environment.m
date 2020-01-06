function [so, poolobj, behavfiles] = trust_setup_environment(so)
%this function sets up the paths for sceptic, VBA, and expected behavior files based on the user's environment
%it also handles setup of the parpool depending on the user

%os = computer;
[~, me] = system('whoami');
me = strtrim(me);
is_alex= strcmp(me,'Alex')==1 || strcmp(me,'dombax')==1;
is_jiazhouchen=strcmp(me,'jiazhouchen')==1;
%is_jiazhou_t=strcmp(me,'jiazhou')==1;
is_mnh=strcmp(me, 'mnh5174') ==1;
is_alison_mbp = strcmp(me, 'alisonmarie526') == 1;
is_alison_mac = strcmp(me, 'ams939') == 1;
is_daniel_mac = strcmp(me,'dpp5430') == 1;
is_polina_mac = strcmp(me, 'polinavanyukov') == 1;

%% set environment and define file locations
if is_alex
    so.vba_path='~/code/VBA-toolbox';
    so.trust_repo='';
    so.output_dir = '~/Box Sync/skinner/projects_analyses/Trust/mfx_analyses';
elseif is_jiazhouchen
    so.vba_path='/Volumes/bek/vbatoolbox';
    so.boxdir = '~/Box';
    so.output_dir = '~/Box/skinner/projects_analyses/Trust/mfx_analyses';
elseif is_mnh
    so.vba_path='~/Documents/MATLAB/VBA-toolbox';
    so.trust_repo='/Users/mnh5174/Data_Analysis/Social_Trust';
    so.data_directory=[so.trust_repo, '/data/', so.dataset, '/mat'];
    so.output_dir=[so.trust_repo, '/output/vba_mfx'];
elseif is_alison_mbp
    so.vba_path='/Users/alisonmarie526/Downloads/VBA-toolbox-master-2';
    so.trust_repo='/Users/alisonmarie526/github_dirs/Social_Trust';
    so.output_dir = [so.trust_repo, '/output/vba_mfx'];
    so.data_directory = [so.trust_repo,'/data/',so.dataset,'/mat'];
    mmrescan_dir = '/Users/alisonmarie526/github_dirs/Social_Trust';
elseif is_alison_mac
    so.vba_path='/Users/ams939/Downloads/VBA-toolbox-master-2';
    so.trust_repo = '~/Box Sync/DEPENd/Projects/Social_Trust';
    so.output_dir = [so.trust_repo, '/output/vba_mfx'];
    mmrescan_dir = '/Users/ams939/Box Sync/DEPENd/Projects/Social_Trust';
elseif is_daniel_mac
    so.vba_path='/Users/dpp5430/Desktop/VBA-toolbox';
    so.trust_repo = '/Users/dpp5430/Desktop/workingST';
    so.output_dir = '/Users/dpp5430/Desktop/output';
    so.data_directory = [so.trust_repo,'/data/',so.dataset,'/mat'];
elseif is_polina_mac
    so.vba_path='/Users/polinavanyukov/Scripts/VBA-toolbox';
    so.boxdir = '/Users/polinavanyukov/Box Sync/';
    so.data_dir = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/';
    so.sbjlist = '/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/trust_scan_masterlist.txt'; %currently not in use
    so.output_dir = ['/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/',so.model,'/',so.subj_subset,'/'];
else    
    so.output_dir = [trust_repo, '/output/vba_out/', so.dataset, '/mfx/', so.model];
    if ~exist(so.output_dir, 'dir'), mkdir(so.output_dir); end
end

addpath(genpath_safe(so.vba_path)); %add VBA functions


if ~isfield(so, 'dataset'), so.dataset='mmrescan'; end %default
if ~isfield(so, 'subj_subset'), so.subj_subset=''; end %fit all



if strcmpi(so.dataset, 'mmrescan')
    %% identify subject behavior files
    demographics_file = readtable([so.trust_repo, '/Rescanning_Project_Scan_Log.csv']);
    %demographics_file = removevars(demographics_file, {'Notes', 'Confidence'});
    matching_file = readtable([so.trust_repo, '/subsetting_pts/mmrescan_matchingresults_rist_jan2019.csv']);
    %matching_file = removevars(matching_file, {'Var1', 'matchid'});
    if strcmpi(so.subj_subset, 'normative')
        %dataset containing all controls
        controls=strcmpi(demographics_file.Dx_Group, 'HC');
        behavfiles=strcat(so.data_directory, '/trust', demographics_file.MR_Center_ID(controls), '.mat');
    elseif strcmpi(so.subj_subset, 'matched')
        behavfiles=strcat(so.data_directory, '/trust', matching_file.MR_Center_ID, '.mat');
    elseif strcmpi(so.subj_subset, 'specc')
        %dataset containing both BPD and matched control participants
        specc=contains(demographics_file.Study, 'SPECC');
        behavfiles=strcat(so.data_directory, '/trust', demographics_file.MR_Center_ID(specc), '.mat');
    else
        fprintf('Fitting all subjects because subj_subset not used/matched\n');
        behavfiles=strcat(so.data_directory, '/trust', demographics_file.MR_Center_ID, '.mat');
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
elseif strcmpi(so.dataset, 'bsocial')
    so.data_directory=fullfile(so.boxdir,'skinner/data/eprime/bpd_trust');
    addpath import;
    IDdirs=dir(so.data_directory);
    IDs={IDdirs.name};
    IDs=IDs(~cellfun(@(x) strcmpi(x,'.') | strcmpi(x,'..'),IDs));
    allfilepath = []; allb=[];
    for id = IDs
        IDtext=sprintf('ID%s', id{:});
        [allb.(IDtext),allfilepath.(IDtext)] = trustbehavior(so.data_directory,id{:},0);
    end
    so.allb=allb;
    so.n_trials=144;
    behavfiles=struct2cell(allfilepath);
    behavfiles=behavfiles(~cellfun(@isempty,behavfiles));
elseif strcmpi(so.dataset,'learn')
    so.data_directory=fullfile(so.boxdir,'skinner/data/eprime/scan_trust');
    addpath import;
    if contains(so.subj_subset,'all')
        IDs=[{'35780'}, {'46069'}, {'202200'}, {'202278'}, {'203264'}, {'203803'}, {'204015'}, {'206270'}, {'207224'}, {'208510'}, {'210290'}, {'210374'}, {'210701'}, {'211038'}, {'211101'}, {'211253'}, {'212019'}, {'212020'}, {'212794'}, {'213163'}, {'213704'}, {'213868'}, {'214082'}, {'215088'}, {'215537'}, {'215622'}, {'216011'}, {'216032'}, {'217008'}, {'217048'}, {'217173'}, {'217293'}, {'217988'}, {'218219'}, {'218572'}, {'218714'}, {'219044'}, {'219089'}, {'219392'}, {'219619'}, {'219658'}, {'219676'}, {'219772'}, {'219809'}, {'219886'}, {'219944'}, {'219949'}, {'219954'}, {'220043'}, {'220051'}, {'220064'}, {'220104'}, {'220152'}, {'220184'}, {'220186'}, {'220244'}, {'220299'}, {'220399'}, {'220451'}, {'220477'}, {'220513'}, {'220523'}, {'220531'}, {'220553'}, {'220566'}, {'220597'}, {'220678'}, {'220684'}, {'220740'}, {'220758'}, {'220789'}, {'220796'}, {'220865'}, {'220889'}, {'220913'}, {'220926'}, {'220927'}, {'220947'}, {'220963'}, {'220976'}, {'220989'}, {'221018'}, {'221036'}, {'221050'}, {'221085'}, {'221099'}, {'221140'}, {'221164'}, {'221181'}, {'221183'}, {'221206'}, {'221212'}, {'221226'}, {'221261'}, {'221273'}, {'221292'}, {'221295'}, {'221298'}, {'881075'}, {'881105'}, {'881106'}, {'881132'}, {'881224'}, {'881230'}];
    elseif contains(so.subj_subset,'HC')
        IDs=[{'202200'},{'210290'},{'210374'},{'211253'},{'212019'},{'212020'},{'212794'},{'213163'},{'213868'},{'214082'},{'215537'},{'216011'},{'216032'},{'219044'},{'219886'},{'219944'},{'220043'},{'220051'},{'220104'},{'220184'},{'220244'},{'220299'},{'220740'},{'221050'},{'221085'},{'221140'},{'221164'},{'221183'},{'221261'}];
    elseif contains(so.subj_subset,'DEP')
        IDs=[{'206270'}, {'207224'}, {'208510'}, {'211038'}, {'213704'}, {'215088'}, {'218714'}, {'219809'}, {'219949'}, {'219954'}, {'220477'}, {'220597'}, {'220789'}, {'220926'}, {'220963'}, {'221212'}, {'221273'}, {'221295'}, {'881105'}, {'881106'}, {'881132'}, {'881224'}, {'881230'}];
    elseif contains(so.subj_subset, 'IDEA')
        IDs=[{'203803'}, {'210701'}, {'217048'}, {'218572'}, {'219658'}, {'219676'}, {'220152'}, {'220553'}, {'220678'}, {'220758'}, {'220796'}, {'220865'}, {'220913'}, {'220947'}, {'220976'}, {'221018'}, {'221036'}, {'221181'}, {'221292'}];
    else %ATT
        IDs=[{'35780'},{'46069'}, {'202278'}, {'203264'}, {'204015'}, {'211101'}, {'215622'}, {'217008'}, {'217173'}, {'217293'}, {'217988'}, {'218219'}, {'219089'}, {'219392'}, {'219619'}, {'219772'}, {'220064'}, {'220186'}, {'220399'}, {'220451'}, {'220513'}, {'220523'}, {'220531'}, {'220566'}, {'220684'},{'220889'}, {'220927'}, {'220989'}, {'221206'}, {'221226'}, {'221298'}, {'881075'}];
    end
    allfilepath = []; allb=[];
    for id = IDs
        IDtext=sprintf('ID%s', id{:});
        [allb.(IDtext),allfilepath.(IDtext)] = trustbehavior(so.data_directory,id{:},0);
    end
    so.allb=allb;
    so.n_trials=192;
    behavfiles=struct2cell(allfilepath);
    behavfiles=behavfiles(~cellfun(@isempty,behavfiles));
end


%% setup parallel parameters
if is_alex
    if strcmp(me,'dombax')==1
        ncpus = 12;
        fprintf('defaulting to 12 cpus on Thorndike \n');
    else
        ncpus=10;
        fprintf('defaulting to 10 cpus on iMac Pro \n');
    end
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
    end
elseif is_jiazhouchen
    ncpus=4;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
    end
else
    ncpus=getenv('matlab_cpus');
    if strcmpi(ncpus, '')
        ncpus=2;
        fprintf('defaulting to 2 cpus because matlab_cpus not set\n');
    else
        ncpus=str2double(ncpus);
    end
    
    poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
end

end
