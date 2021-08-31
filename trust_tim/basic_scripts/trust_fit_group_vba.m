%loads in subjects' data and fits SCEPTIC models using VBA;
close all;
clear;

vo=[]; %VBA options object
vo.trustees = {'neutral', 'bad', 'good'}; %used for setting up dummy code scheme for kappaT. All subjects need the same trustee set
vo.fit_groupfixed=true; %whether to refit model at group posterior estimates
vo.save_all = true;

%model nomenclature:
% sharing bias parameters:
%   'nK' - no Kappa
%   'sK' - subject only kappa
%   'stK' - subject and trustee kappa

%left side of pipe is model (subfolder name), right side of pipe is the variant, including any additional
%parameterization details like subject-specific kappa.

%vo.model = {'policy|stK_mfx', 'twoQs|Actual_Outcome_stK_ffx', 'rational|rational_coopbonus_sK_mfx'};

%set models
%vo.models_to_run = {'twoQs|Actual_Outcome_stK_ffx', 'rational|rational_stK_ffx'}; %basal rational model
%vo.models_to_run = {'policy|stK_mfx'}
vo = trust_validate_options(vo); %initialize and validate suuvid fitting settings

%vo.subj_subset = 'normative';  %which subset of files to fit
%vo.subj_subset = 'matched';
%vo.subj_subset = '';

if ~strcmpi(vo.subj_subset, '')
  dataset_name=[vo.dataset, '_', vo.subj_subset];
else
  dataset_name=vo.dataset;
end 

[vo, poolobj, behavfiles] = trust_setup_environment(vo);

vba_working_dir = fileparts(mfilename('fullpath')); %for paths relative to the vba repo
models_dir = [ vba_working_dir, filesep, 'models' ]; %contains compartmentalized code for each model
addpath(models_dir); %shared functions across models can live in models root

% test selective wide maintenance/stick PE and uniform stick PE
%vo.models_to_run={'policy_nomult_stK', 'policy_nomult_kappaSonly'};
%vo.models_to_run = {'policy_mult_fixed4_stK', 'policy_nomult_stK_ffx', 'policy_nomult_kappaSonly_ffx', 'policy_nomult_kappaSonly', 'policy_mult_stK'};
%vo.models_to_run = {'policy_nomult_stK_ffx', 'policy_nomult_stK'};
%vo.models_to_run = {'mmrescan_countertrustee_nomult_stK_ffx', 'mmrescan_countertrustee_nomult_stK', 'mmrescan_regret_nomult_stK', 'mmrescan_regret_nomult_stK_ffx', 'mmrescan_null_nomult_stK_ffx', 'mmrescan_null_nomult_stK'};
%vo.models_to_run = {'mmrescan_prs_nomult_stK'};
%vo.models_to_run = {'policy|policy_stK_mfx_coopbonus'};
%vo.models_to_run = {'rational|pavvalue_stK_mfx'};

for i = 1:length(vo.models_to_run)
    vo.model = vo.models_to_run{i}; %set model string for this model
    ffx = contains(vo.model, 'ffx'); %mfx vs. ffx
    if ffx, vo.fitting_procedure='ffx'; else, vo.fitting_procedure='mfx'; end
    fprintf('fitting model: %s\n', vo.model);
    if contains(vo.model, '|')
        ss = split(vo.model, '|');
        vo.model = ss{1};
        vo.variant = ss{2}; %variant of this model (all parameters the same, just the equations differ)
    else
        vo.variant = '';
    end
    
    modelinfo = strrep(vo.models_to_run{i}, '|', '_'); %for naming output
    this_odir = [vo.output_dir, filesep, vo.fitting_procedure, filesep, dataset_name, filesep, vo.model]; %output directory for this model
    if ~exist(this_odir, 'dir'), mkdir(this_odir); end
    
    vo.graphics = 0; %don't display fitting interactively
    vo.multisession = 0; %not currently supported; kappa trustee handled by dummy coding
    vo.fixed=0;
    vo.censor=0;
    
    %setup this model by adding its parts to the path
    mdir = [models_dir, filesep, vo.model];
    if ~exist(mdir, 'dir'), error('cannot locate model directory: %s', mdir); else, addpath(mdir); end
    vo = m_setup(vo);
    
    %vo.params = 'fixedparams';
    vo.params= 'freed';
    
    %extract IDs for record keeping
    % [~,fnames]=cellfun(@fileparts, behavfiles, 'UniformOutput', false);
    % ids=cellfun(@(x) char(regexp(x,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match')), fnames, 'UniformOutput', false);
    
    %puts a nested cell in each element
    %ids=regexp(fnames,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)
    
    %number of subjects to fit
    
    ns = length(behavfiles);
    y_all = cell(ns, 1);
    u_all = cell(ns, 1);
    options_all = cell(ns, 1);
    hasComputer_all = cell(ns, 1);
    noComputer_all = cell(ns, 1);
    n_t = NaN(1,ns);
    ids = cell(1, ns);
    for sub = 1:ns
        [~,id,~] = fileparts(behavfiles{sub});
        id=erase(id, 'trust'); %remove 'trust' from id string
        ids{sub} = id;
        fprintf('Loading subject %d id: %s \n', sub, id);
        [y, u, noresponse, hasComputer] = trust_get_data(behavfiles{sub}, vo);
        vo.n_t = length(y);
        [options, dim] = mmrescan_get_vba_options(vo);
        n_t(sub) = dim.n_t; % allow for variation in number of trials across subjects
        %options.skipf = noresponse | options.skipf;
        % populate data structures for VBA_MFX
        hasComputer_all{sub} = hasComputer;
        noComputer_all{sub} = ~hasComputer;
        y_all{sub} = y;
        u_all{sub} = u;
        options_all{sub} = options;
    end
    
    %     if (PConly)
    %         hasComputer_all=cell2mat(hasComputer_all);
    %         nomulttext=strcat(nomulttext,'_PCOnly');
    %         ids = ids(hasComputer_all);
    %         n_t = n_t(hasComputer_all);
    %         y_all = y_all(hasComputer_all);
    %         u_all = u_all(hasComputer_all);
    %         options_all = options_all(hasComputer_all);
    %         ns = length(y_all);
    %     end
    %
    %     if (NoCompOnly)
    %         noComputer_all=cell2mat(noComputer_all);
    %         nomulttext=strcat(nomulttext,'_noPCOnly');
    %         ids = ids(noComputer_all);
    %         n_t = n_t(noComputer_all);
    %         y_all = y_all(noComputer_all);
    %         u_all = u_all(noComputer_all);
    %         options_all = options_all(noComputer_all);
    %         ns = length(y_all);
    %     end
    
    % add the n_t vector to dim before passing to MFX
    dim.n_t=n_t;
    dim.n_s=length(y_all);
    
    %options for MFX
    options_group.TolFun=1e-2;
    options_group.MaxIter=50;
    options_group.DisplayWin=1;
    options_group.verbose=1;
    priors_group = options_all{1}.priors; %use first subject as representative of group priors on phi, theta, and X0
    
    %priors_group.muPhi = zeros(dim.n_phi,1); %temperature -- exp(phi(1))
    %priors_group.SigmaPhi = 1e1*eye(dim.n_phi); %variance on temperature (before exponential transform)
    %priors_group.muTheta = zeros(dim.n_theta,1); %learning rate (alpha), selective maintenance (gamma) -- before logistic transform
    %priors_group.SigmaTheta =  1e1*eye(dim.n_theta); %variance of 10 on alpha and gamma
    %priors_group.muX0 = zeros(vo.hidden_states,1); %have PE and decay as tag-along states
    %priors_group.SigmaX0 = zeros(vo.hidden_states, vo.hidden_states); %have PE and decay as tag-along states
    
    priors_group.a_vX0 = Inf([1, vo.hidden_states]); %use infinite precision prior on gamma for X0 to treat as fixed (a = Inf; b = 0)
    priors_group.b_vX0 = zeros([1, vo.hidden_states]);
    dim_share = dim;
    
    %% core parameter fitting loop, either using FFX single-subject fitting or MFX fitting
    if ffx
        parfor s = 1:ns
            dim_s = dim_share;
            dim_s.n_t = n_t(s);
            fprintf('Running %d\n',s);
            [p_sub_tmp, o_sub_tmp] = VBA_NLStateSpaceModel(y_all{s}, u_all{s}, vo.evo_fname, vo.obs_fname, dim_s, options_all{s});
            
            p_sub{s} = p_sub_tmp;
            o_sub{s} = o_sub_tmp;
        end
        
        for s = 1:ns
            o_sub{s}.options.inF.id = ids{ns};
        end
    else
        
        [p_sub, o_sub, p_group, o_group] = VBA_MFX_parallel(y_all, u_all, vo.evo_fname, vo.obs_fname, dim, options_all, priors_group, options_group);
    end
    
    %populate subject ids into output.options.inF structure since these are added to subject statistics below
    for s=1:length(o_sub)
        o_sub{s}.options.inF.id = ids{s};
        o_sub{s}.options.inG.id = ids{s};
        
        %populate ffx parameters from o_group structure to p_sub structure for extraction
        if ~ffx
            p_sub{s}.ffx = o_group.initVBA.p_sub{s};
            p_sub{s}.ffx = rmfield(p_sub{s}.ffx, {'SigmaX', 'iQx', 'iQy'}); %large matrices not needed for anything (they use lots of disk space)
        end
    end
    
    %% save model outputs
    if vo.save_all
        %too huge to save into one .mat file
        save([this_odir, '/', dataset_name, sprintf('_%s_vba_results_psub.mat', modelinfo)], 'p_sub','-v7.3');
        save([this_odir, '/', dataset_name, sprintf('_%s_vba_results_osub.mat', modelinfo)], 'o_sub','-v7.3');
        
        if ~ffx %save group posteriors for mfx fitting
            save([this_odir, '/', dataset_name, sprintf('_%s_vba_results_pgroup.mat', modelinfo)], 'p_group','-v7.3');
            save([this_odir, '/', dataset_name, sprintf('_%s_vba_results_ogroup.mat', modelinfo)], 'o_group','-v7.3');
        end
    end
    
    %save log evidence vector for model
    if ~ffx
        F = o_group.within_fit.F;
    else
        F = zeros(ns, 1);
        for xx = 1:ns
            F(xx) = o_sub{xx}.F;
        end
    end
    
    save([this_odir, '/', dataset_name, sprintf('_%s_vba_F.mat', modelinfo)], 'F');
    
    %create a structure with just the barebones useful parameters from subjects
    s_all = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub, o_sub, 'UniformOutput', false);
    
    %output global statistics and trial-level outputs as CSVs    
    [group_global, group_trial_level] = extract_group_statistics(s_all, ...
        sprintf('%s/%s_%s_global_statistics.csv', this_odir, dataset_name, modelinfo), ...
        sprintf('%s/%s_%s_trial_outputs.csv', this_odir, dataset_name, modelinfo));
       
    %save group statistics and model-fitting inputs
    save(sprintf('%s/%s_%s_group_fits', this_odir, dataset_name, modelinfo), 'ids', 'vo', 's_all', 'group_global', 'group_trial_level');
    save(sprintf('%s/%s_%s_vba_results_settings.mat', this_odir, dataset_name, modelinfo), 'priors_group', 'options_group', 'y_all', 'u_all', 'ids', '-v7.3');

    %% refit subject data at group means from MFX -- 'group-fixed' approach
    if ~ffx
        phi_fix = p_group.muPhi;
        theta_fix = p_group.muTheta;
    else
        mm = cellfun(@(p) p.muPhi, p_sub, 'UniformOutput', false);
        mm = cat(2, mm{:}); %concatenate vectors into a phi x subjects matrix
        phi_fix = mean(mm, 2);
        
        mm = cellfun(@(p) p.muTheta, p_sub, 'UniformOutput', false);
        mm = cat(2, mm{:}); %concatenate vectors into a phi x subjects matrix
        theta_fix = mean(mm, 2);
    end
    
    parfor sub = 1:ns
        fprintf('Running %d through group-fixed fitting\n',sub);
        dim_s = dim;
        dim_s.n_t = n_t(sub);
        
        [p_sub_tmp, o_sub_tmp]=VBA_refit_fixed(y_all{sub}, u_all{sub}, vo.evo_fname, vo.obs_fname, ...
            dim_s, options_all{sub}, phi_fix, theta_fix);
        
        o_sub_tmp.options.id=ids{sub};
        o_sub_tmp.options.inF.id=ids{sub};
        p_sub_fixed{sub} = p_sub_tmp;
        o_sub_fixed{sub} = o_sub_tmp;
    end
    
    s_all_fixed = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub_fixed, o_sub_fixed, 'UniformOutput', false);
    
    %output group-fixed results
    [group_global_fixed, group_trial_level_fixed] = extract_group_statistics(s_all_fixed, ...
        sprintf('%s/%s_%s_global_statistics_groupfixed.csv', this_odir, dataset_name, modelinfo), ...
        sprintf('%s/%s_%s_trial_outputs_groupfixed.csv', this_odir, dataset_name, modelinfo));
    
    save(sprintf('%s/%s_%s_group_fits_groupfixed', this_odir, dataset_name, modelinfo), 'ids', 'vo', 's_all', 'group_global_fixed', 'group_trial_level_fixed');
    
    rmpath(mdir); %cleanup model path
    
end

delete(gcp('nocreate')); %stop any pool that is running
