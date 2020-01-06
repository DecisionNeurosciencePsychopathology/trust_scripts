%loads in subjects' data and fits SCEPTIC models using VBA;
close all;
clear;
%shut down parrallel pool, if one is still running
delete(gcp('nocreate'));

%curpath = fileparts(mfilename('fullpath'));

so=[]; %setup object
so.model='mixe2'; %possible values: policy, mixed, rcounter, rcounte2, rcounte3, mixe2
so.dataset='learn';   %dataset id participants are listed in trust_setup_environment.m
so.subj_subset='all'; %possible values: all, HC, DEP, IDEA, ATT


try
        isqsub=getenv('isqsub');
catch
        isqsub = 0;
end
disp(isqsub);
%if the program is running on the cluster, assign certain variables with the arguments passed in to the script that's being run. Otherwise,these variables can be set manually.
if isqsub
    disp("in");
    so.model=getenv('model');
    so.dataset=getenv('dataset');
    so.subj_subset=getenv('subset');
    modelstorun={getenv('modelstorun')};
    params=getenv('params');
else
    %so.dataset='learn';
    %modelstorun={'policy_nomult_kappaST_mfx'};
    %modelstorun = {'mixed_nomult_kappaST_mfx'};
    %modelstorun = {'rcounte2_nomult_kappaST_mfx'};
    %modelstorun = {'rcounte3_nomult_kappaST_mfx'};
    modelstorun = {'mixe2_nomult_kappaST_mfx'};
    params = 0;
end
disp(so.model)
disp(so.dataset)
disp(so.subj_subset)
disp(modelstorun)
%disp(params)

%so.subj_subset='normative';%which subset of files to fit
%so.subj_subset = 'matched';

[so, poolobj, behavfiles] = trust_setup_environment(so);

% test selective wide maintenance/stick PE and uniform stick PE
%modelstorun={'policy_nomult_kappaST', 'policy_nomult_kappaSonly'};
%modelstorun = {'policy_mult_fixed4_kappaST', 'policy_nomult_kappaST_ffx', 'policy_nomult_kappaSonly_ffx', 'policy_nomult_kappaSonly', 'policy_mult_kappaST'};
save_all = 1;
%modelstorun = {'policy_nomult_kappaST_ffx', 'policy_nomult_kappaST'};
%modelstorun = {'mmrescan_countertrustee_nomult_kappaST_ffx', 'mmrescan_countertrustee_nomult_kappaST', 'mmrescan_regret_nomult_kappaST', 'mmrescan_regret_nomult_kappaST_ffx', 'mmrescan_null_nomult_kappaST_ffx', 'mmrescan_null_nomult_kappaST'};
%modelstorun = {'mmrescan_prs_nomult_kappaST'};


for i = 1:length(modelstorun)
    
    model=modelstorun{i};
    policy=contains(model,'policy');
    mixed = contains(model,'mixed');
    rcounter = contains(model,'rcounter');
    rcounter2 = contains(model,'rcounte2');
    rcounter3 = contains(model,'rcounte3');
    mixed2 = contains(model,'mixe2');
    nomult=contains(model,'nomult');
    fixed1 = contains(model, 'fixed1');
    fixed2 = contains(model, 'fixed2');
    fixed3 = contains (model, 'fixed3');
    fixed4 = contains(model, 'fixed4');
    kappaST=contains(model,'kappaST');
    regret = contains(model, 'regret');
    countertrustee = contains(model, 'countertrustee');
    null = contains(model, 'null');
    kappaSonly = contains(model, 'kappaSonly');
    ED1 = contains(model, 'ED1');
    ED2 = contains(model, 'ED2');
    softmax = contains(model, 'softmax');
    censor = contains(model, 'censor');
    ffx = contains(model, 'ffx');
    mmrescan = contains(model, 'mmrescan');
    ptrustee = contains(model, 'ptrustee');
    utilitygain = contains(model, 'utilitygain');
    prs = contains(model, 'prs');
    PConly = contains(model, 'PConly');
    NoCompOnly = contains(model, 'NoCompOnly');
    fourTrustee = contains(model, 'fourTrustee');
    so.fourTrustee = contains(model, 'fourTrustee');
    %note that this function looks for 'mmrescan_dataset' and 'mmrescan_model'
    %as environment variables so that this script can be scaled easily for batch processing
    
    %% evolution functions
    if (policy)
        so.counter = 2;
        so.model = 'policy';
        infotext=string('policy');
    elseif (mixed)
        so.counter = 2;
        so.model = 'mixed';
        infotext=string('mixed');
    elseif (mixed2)
        so.counter = 2;
        so.model = 'mixed2';
        infotext=string('mixed2');
    elseif (rcounter)
        so.counter = 2;
        so.model = 'rcounter';
        infotext=string('rcounter');
    elseif (rcounter2)
        so.counter = 2;
        so.model = 'rcounter2';
        infotext=string('rcounter2');
    elseif (rcounter3)
        so.counter = 2;
        so.model = 'rcounter3';
        infotext=string('rcounter3');
    elseif (ptrustee)
        so.counter = 2;
        so.model = 'ptrustee';
        infotext = string('ptrustee');
    elseif (prs)
        so.counter = 2;
        so.model = 'prs';
        infotext = string('prs');
    elseif (regret)
        so.model = 'regret';
        so.counter = 1;
        infotext=string('regret');
    elseif (countertrustee)
        so.counter = 3;
        so.model = 'countertrustee';
        infotext = string('countertrustee');
    elseif (utilitygain)
        so.counter = 2;
        so.model = 'utilitygain';
        so.utility_goodbad = contains(model, 'goodbad');
        infotext = string('utilitygain');
    elseif (null)
        so.model = 'null';
        so.counter = 0;
        infotext = string('null');
    end
    
    % multisession
    if (nomult)
        so.multisession = 0;
        so.fixed = 0;
        nomulttext=string('no_multisession');
    elseif (fixed1)
        so.multisession = 1;
        so.fixed = 1;
        mult = 1;
        nomulttext = string('multisession_fixed1');
    elseif (fixed2)
        so.multisession = 1;
        so.fixed = 2;
        mult = 1;
        nomulttext = string('multisession_fixed2');
    elseif (fixed3)
        so.multisession = 1;
        so.fixed = 3;
        nomulttext = string('multisession_fixed3');
    elseif (fixed4)
        so.multisession = 1;
        so.fixed = 4;
        nomulttext = string('multisession_fixed4');
    end
    
    if (~nomult)
        so.n_trustees = 3;
    end
    % PARAMETERIZATION; factorize decay/learning rate (1), don't (0)
    if (kappaST) && (nomult)
        if strcmpi(so.dataset,'learn') || fourTrustee
            so.sigmakappa=5;
            so.fourTrustee = 1;
            factortext=string('FourTrusteekappaST');
        else
            so.sigmakappa=4;
            factortext=string('kappaST');
        end
    elseif (kappaST) && (~nomult)
            so.sigmakappa = 2;
            factortext=string('kappaST_mult');
    elseif (kappaSonly) && (nomult)
        so.sigmakappa = 5;
        so.fourTrustee = 0;
        factortext=string('kappaSonly_nomult');
    elseif (ED1) || (ED2)
        so.sigmakappa = 1;
        factortext=string('kappaSonly_mult');
    elseif (softmax)
        so.sigmakappa = 0;
        factortext = string('noKappa');
    end
    
    if (censor)
        so.censor = 1;
    else
        so.censor = 0;
    end
    
    if isqsub
	params = getenv(params);
    end
    
    if (params)
        so.params = 'fixedparams';
    else
        so.params = 'freed';
    end
    
    so = mmrescan_validate_options(so); %initialize and validate fitting settings
    
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
    n_t=NaN(1,ns);
    ids = cell(1, ns);
    for sub = 1:ns
        [~,id,~] = fileparts(behavfiles{sub});
        id=erase(id, 'trust'); %remove 'trust' from id string
        ids{sub} = id;
        fprintf('Loading subject %d id: %s \n', sub, id);
        [y, u, noresponse,hasComputer] = mmrescan_get_data(behavfiles{sub}, so);
        so.n_t = length(y);
        [options, dim] = mmrescan_get_vba_options(so);
        n_t(sub) = dim.n_t; % allow for variation in number of trials across subjects
        %options.skipf = noresponse | options.skipf;
        % populate data structures for VBA_MFX
        hasComputer_all{sub} = hasComputer;
        noComputer_all{sub} = ~hasComputer;
        y_all{sub} = y;
        u_all{sub} = u;
        options_all{sub} = options;
    end
    
    if(PConly)
        hasComputer_all=cell2mat(hasComputer_all);
        nomulttext=strcat(nomulttext,'_PCOnly');
        ids = ids(hasComputer_all);
        n_t = n_t(hasComputer_all);
        y_all = y_all(hasComputer_all);
        u_all = u_all(hasComputer_all);
        options_all = options_all(hasComputer_all);
        ns = length(y_all);
    end
    
    if(NoCompOnly)
        noComputer_all=cell2mat(noComputer_all);
        nomulttext=strcat(nomulttext,'_noPCOnly');
        ids = ids(noComputer_all);
        n_t = n_t(noComputer_all);
        y_all = y_all(noComputer_all);
        u_all = u_all(noComputer_all);
        options_all = options_all(noComputer_all);
        ns = length(y_all);
    end
    
    % add the n_t vector to dim before passing to MFX
    dim.n_t=n_t;
    dim.n_s=length(y_all);
    
    %options for MFX
    options_group.TolFun=1e-2;
    options_group.MaxIter=50;
    options_group.DisplayWin=1;
    options_group.verbose=1;
    
    priors_group.muPhi = zeros(dim.n_phi,1); %temperature -- exp(phi(1))
    priors_group.SigmaPhi = 1e1*eye(dim.n_phi); %variance on temperature (before exponential transform)
    priors_group.muTheta = zeros(dim.n_theta,1); %learning rate (alpha), selective maintenance (gamma) -- before logistic transform
    priors_group.SigmaTheta =  1e1*eye(dim.n_theta); %variance of 10 on alpha and gamma
    priors_group.muX0 = zeros(so.hidden_states,1); %have PE and decay as tag-along states
    priors_group.SigmaX0 = zeros(so.hidden_states, so.hidden_states); %have PE and decay as tag-along states
    priors_group.a_vX0 = Inf([1, so.hidden_states]); %use infinite precision prior on gamma for X0 to treat as fixed (a = Inf; b = 0)
    priors_group.b_vX0 = zeros([1, so.hidden_states]);
    dim_share = dim;
    if (~nomult) || (ffx)
        for s= 1:ns
            dim_s = dim_share;
            dim_s.n_t = n_t(s);
             fprintf('Running %d\n',s);
            [p_sub_tmp, o_sub_tmp] = VBA_NLStateSpaceModel(y_all{s}, u_all{s}, so.evo_fname, so.obs_fname, dim_s, options_all{s});
            
            p_sub{s} = p_sub_tmp;
            o_sub{s} = o_sub_tmp;
        end
        for s = 1:ns
            o_sub{s}.options.inF.id = ids{ns};
        end
    else
         [p_sub, o_sub, p_group, o_group] = VBA_MFX_parallel(y_all, u_all, so.evo_fname, so.obs_fname, dim, options_all, priors_group, options_group);
    end
    delete(poolobj);
    
    %populate subject ids into output.options.inF structure since these are added to subject statistics below
    if (nomult) && (~ffx)
        for s=1:length(o_sub)
            o_sub{s}.options.inF.id = ids{s};
            o_sub{s}.options.inG.id = ids{s};
            
            %populate ffx parameters from o_group structure to p_sub structure for extraction
            p_sub{s}.ffx = o_group.initVBA.p_sub{s};
            p_sub{s}.ffx = rmfield(p_sub{s}.ffx, {'SigmaX', 'iQx', 'iQy'}); %large matrices not needed for anything (they use lots of disk space)
        end
        
        if save_all            
            %sprintf('_%s_%s_%s_vba_mfx_results_settings.mat', infotext, tempgentext, factortext)
            
            %too huge to save into one .mat file
            save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_mfx_results_psub.mat', infotext, nomulttext, factortext)], 'p_sub');
            save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_mfx_results_pgroup.mat', infotext, nomulttext, factortext)], 'p_group');
            save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_mfx_results_ogroup.mat', infotext, nomulttext, factortext)], 'o_group');
            save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_mfx_results_osub.mat', infotext, nomulttext, factortext)], 'o_sub');
            save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_mfx_results_settings.mat', infotext, nomulttext, factortext)], 'priors_group', 'options_group', 'y_all', 'u_all', 'ids');
        end
        F = o_group.within_fit.F;
        % just the L
        save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_mfx_L.mat', infotext, nomulttext, factortext)], 'F');
      
        %create a structure with just the barebones useful parameters from subjects
        s_all = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub, o_sub, 'UniformOutput', false);
      
        %output MFX results
        infotext='MFXPara';
        [group_global, group_trial_level] = extract_group_statistics(s_all, ...
            sprintf('%s/%s_%s_%s_%s_%s_mfx_global_statistics.csv', so.output_dir, so.dataset, so.model, infotext, nomulttext, factortext), ...
            sprintf('%s/%s_%s_%s_%s_%s_mfx_trial_outputs_by_timestep.csv', so.output_dir, so.dataset, so.model, infotext, nomulttext, factortext));
        
        %fixed parameters (GetPosterior_FixPara to run everything w/
        %infinite precision priors
        %group median parameters

        for sub= 1:ns
            fprintf('Running %d\n',sub);
            options_all{sub}.dim = dim;
            options_all{sub}.dim.n_t = options_all{sub}.dim.n_t(sub);
            options_all{sub}.f_fname = so.evo_fname;
            options_all{sub}.g_fname = so.obs_fname;
            options_all{sub}.g_fname = so.obs_fname;
            [p_sub_tmp,o_sub_tmp]=GetPosterior_FIXPara(u_all{sub},y_all{sub},options_all{sub},p_group.muPhi,p_group.muTheta);
            o_sub_tmp.options.id=ids{sub};
            o_sub_tmp.options.inF.id=ids{sub};
            p_sub{sub} = p_sub_tmp;
            o_sub{sub} = o_sub_tmp;
        end

        s_all = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub, o_sub, 'UniformOutput', false);
        
        %output MFX results
        infotext='FixedPara';
        [group_global, group_trial_level] = extract_group_statistics(s_all, ...
            sprintf('%s/%s_%s_%s_%s_%s_mfx_global_statistics.csv', so.output_dir, so.dataset, so.model, infotext, nomulttext, factortext), ...
            sprintf('%s/%s_%s_%s_%s_%s_mfx_trial_outputs_by_timestep.csv', so.output_dir, so.dataset, so.model, infotext, nomulttext, factortext));
            
        
        
        %save the basis as a csv file
        if any(strcmpi(so.model, {'policy', 'countertrustee', 'null', 'regret', 'ptrustee','utilitygain' ,'prs','mixed','rcounter','rcounter2', 'rcounter3','mixed2'}))
            so.ntimesteps = 1;
        end
        vnames = cellfun(@(x) strcat('Time_', num2str(x)), num2cell(1:so.ntimesteps), 'UniformOutput', false);
        %basis_mat = array2table(o_sub{1}.options.inF.gaussmat, 'VariableNames', vnames);
        
         %save group outputs
         save(sprintf('%s/group_fits_%s_%s_%s_%s_%s', so.output_dir, so.model, so.dataset,infotext, nomulttext, factortext), 'ids', 'so', 's_all', 'group_global', 'group_trial_level');
          
       
        
    else
        s_all = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub, o_sub, 'UniformOutput', false);
        
        %output MFX results
 
        [group_global, group_trial_level] = extract_group_statistics(s_all, ...
            sprintf('%s/%s_%s_%s_%s_%s_mfx_global_statistics.csv', so.output_dir, so.dataset, so.model, infotext, nomulttext, factortext), ...
            sprintf('%s/%s_%s_%s_%s_%s_mfx_trial_outputs_by_timestep.csv', so.output_dir, so.dataset, so.model, infotext, nomulttext, factortext));
        
        save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_ffx_results_psub.mat', infotext, nomulttext, factortext)], 'p_sub');
        save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_ffx_results_osub.mat', infotext, nomulttext, factortext)], 'o_sub');
        negF = zeros(ns, 1);
        for xx = 1:ns
            negF(xx) = o_sub{xx}.F;
        end
        save([so.output_dir, '/', so.dataset, '_', so.model, sprintf('_%s_%s_%s_vba_ffx_L.mat', infotext, nomulttext, factortext)], 'negF');
        
    end
    
    
end
