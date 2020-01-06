clear variables;
close all;

%% setting data location + masterlist
scanner_subjs = 1; %Set to 1 is scanning subjects processed, 0 is behavioral
clinical = 0; %for clinical datasets
hallquist_override = 0; %If using Hallquist replication dataset
save_str = date; %default

os = computer;
if strcmp(os(1:end-2),'PCWIN')
    path_to_trust_ids = 'C:/kod/trust_basic_scripts\';
    if scanner_subjs
        datalocation = glob('E:/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        masterlist = load([path_to_trust_ids 'trust_ids']); %Need a way to automatically update these text files!
    else
        datalocation = glob('E:/Box Sync/Project Trust Game/data/processed/beha_behavior/');
        masterlist = load([path_to_trust_ids 'trust_ids_behav']); %Need a way to automatically update these text files!
        save_str = 'final_behavior/';
    end
else
    [~, me] = system('whoami');
    me = strtrim(me);    
    if strcmp(me,'polinavanyukov')==1
        if clinical == 0
            path_to_trust_ids = '/Users/polinavanyukov/Scripts/Trust/trust_rl_VBA/';
            if scanner_subjs == 0
               datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/');
               masterlist = load([path_to_trust_ids 'trust_ids_behav']); %behavioral master list
            else
               datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
               masterlist = load([path_to_trust_ids 'trust_ids']);    %scan master list
            end
        else
            path_to_trust_ids = '/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/';
            if scanner_subjs == 0
               datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/');
               masterlist = tdfread([path_to_trust_ids 'trust_beha_masterlist.txt']); %behavioral master list
               masterlist.trust_ids = masterlist.subjectID;
               save_str = 'clinical_behaviorSEP04/';
            else
               datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
               masterlist = tdfread([path_to_trust_ids 'trust_scan_masterlist.txt']);    %scan master list
               masterlist.trust_ids = masterlist.subjectID;
               save_str = 'clinical_scanSEP04/';
            end
        end
    else
        datalocation = glob('?');
    end
end

if hallquist_override
    %path to trust ids is the same currently
    %datalocation = glob('C:/Users/wilsonj3/Desktop/hallquist_trust/new_data_1_13/');
    datalocation = glob('/Users/polinavanyukov/Documents/scripts/Trust/hallquist_01132017/');
    %masterlist = load([path_to_trust_ids 'hallquist_trust_ids']); %Need a way to automatically update these text files!
    %masterlist = load([path_to_trust_ids 'hallquist_trust_HCs']); 
    masterlist = load([path_to_trust_ids 'hallquist_trust_HC_03012017']); 
%    save_str = 'new_';
end
%% chose models to fit, defaults to Q-learning
%modelnames = {'ushifted_trust_Qlearning'};

%% set parameters
censor = 0;                     %optional censoring of computer block
sigma_kappa = 2;                %kappa (or action bias) parameter; kappa = 0: no bias; kappa = 1: single bias parameter; kappa = 2: subject + trustee-specific bias parameters
counter = 2;                    %counter = 0: no counterfactual feedback; counter = 1: regret; counter = 2: subject/policy; 3: trustee; 4: other type of regret (reviewers); else counter = 5: disappointment (reviewers); counter = 6: SVM model
multisession = 1;               %0 = modelling all runs as a single block; 1 = modelling runs separately
fixed = 3;                      % set multisession priors (see description below)
%1: X0 is free; subject-wise kappa only; 
%2: X0 is fixed; trustee-wise kappa is free; 
%3: X0 fixed/subject-wise kappa fixed/condition-level kappa is free; 
%4: X0 is fixed; subject-wise kappa only; learning rate (alpha) is free;
%5: X0 is fixed; subject-wise kappa only; learning-rate is fixed; beta is
%free.
%6: All parameters are fixed to a particular constant value (for generating
%value and PE regressors with group median parameter values).


err_counter = 1;

%% wrapper loop
for index=1:length(masterlist.trust_ids)
    id = masterlist.trust_ids(index);
    %% Qlearning with multisession kappa parameters
    if counter ~= 6 %not SVM model
        try
            [posterior, out] = trust_Qlearning_ushifted_v2(datalocation{:}, id, counter, multisession,...
                fixed, sigma_kappa, censor, save_str);
            L(index) = out.F;
        catch ME
            [posterior, out] = NaN;
            error_id(err_counter) = id;
            err_counter = err_counter+1;
        end
    else
        % add on for SVM (need to re-write this clunky bit)
        if (id ~= 215129 && id ~= 220116 && id ~= 219944) %these participants do not have survey results
            [posterior, out] = trust_Qlearning_ushifted_SVM(datalocation{:}, id, counter, multisession,...
                fixed, sigma_kappa ,save_str,1);
             L(index) = out.F;
        end
    end
end

L_name = sprintf('L_cntr%d_mltrun%d_fixed%d_kappa%d_censor%d',counter, multisession,fixed, sigma_kappa, censor);
save(char(L_name), 'L');