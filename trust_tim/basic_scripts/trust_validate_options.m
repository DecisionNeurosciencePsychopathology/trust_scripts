function [vo] = trust_validate_options(vo)
% nbasis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
% fixed_params_across_runs -- self-explanatory
% fit_propspread -- makes temporal generalization within the eligibility trace a free parameter
% ntimesteps:      number of time bins
% u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum

if nargin < 1, vo=[]; end

%if user specifies a dataset upstream, don't check environment variable
if isfield(vo, 'dataset')
    fprintf('Using user-specified dataset: %s.\n   Will not check ''dataset'' environment variable', vo.dataset);
else
    vo.dataset=getenv('dataset');
    if strcmpi(vo.dataset, '')
        %vo.dataset='mmrescan';
         vo.dataset='bsocial';
    end
end

%if user specifies a subset upstream, don't check environment variable
if isfield(vo, 'subj_subset')
    fprintf('Using user-specified subset: %s.\n   Will not check ''subj_subset'' environment variable', vo.subj_subset);
else
    vo.subj_subset=getenv('subj_subset');
    if strcmpi(vo.subj_subset, '')
        vo.subj_subset=''; %defaults to all subjects
    end
end

%same principle for model variant
if isfield(vo, 'models_to_run')
    fprintf('Using user-specified model string: %s.\n   Will not check ''models_to_run'' environment variable', strjoin(vo.models_to_run, ', '));
else
    vo.models_to_run=getenv('models_to_run');
    if strcmpi(vo.models_to_run, '')
        vo.models_to_run = {'policy|stK_mfx'}; % policy model with trustee kappas and MFX estimation from Vanyukov 2019
    else
        vo.models_to_run = strsplit(vo.models_to_run, ',');
    end
end

if ~isfield(vo, 'multisession'), vo.multisession=0; end
if ~isfield(vo, 'graphics'), vo.graphics=0; end
if ~isfield(vo, 'n_trials'), vo.n_trials=144; end
if ~isfield(vo, 'saveresults'), vo.saveresults=0; end %used in some places to denote whether to save fitted outputs

vo.trials_per_run = 48;

%vo.n_t = vo.n_trials;
% vo.n_hidden_states = vo.hidden_states;

%old vestiges of model setup for variants -- now abstracted to models folder

%hidden states field used to index hidden state vector inside evolution and observation functions
% if strcmpi(vo.model,'rationalsd') | strcmpi(vo.model,'rationalpol') | strcmpi(vo.model,'rationalreg') | strcmpi(vo.model,'rationalremod') | strcmpi(vo.model,'rationalcounter') | strcmpi(vo.model,'rationalfullpol') | strcmpi(vo.model, 'ratcoop') | strcmpi(vo.model, 'ratrightness')
%     vo.evo_fname = @f_trust_Qlearn_utilitygain;
%     vo.hidden_states=2;
%     vo.state_names = {'P_share', 'P_PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model,'utilitygain')
%     vo.evo_fname = @f_trust_Qlearn_utilitygain;
%     vo.hidden_states=2;
%     if vo.utility_goodbad
%         vo.n_theta=2;
%         vo.theta_names={'alpha_good','alpha_bad'};
%     else
%         vo.n_theta=1;
%         vo.theta_names={'alpha'};
%     end
%     vo.state_names = {'P_share', 'P_PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model,'ratrightness')
%     vo.n_theta = 1;
%     vo.theta_names = {'alpha'};
% elseif strcmpi(vo.model,'policy')
%     vo.evo_fname = @f_trust_Qlearn_policy;
%     vo.hidden_states=2; %track basis weights (value), as well as PE as tag-along state
%     vo.n_theta=1; %learning rate
%     vo.theta_names={'alpha'};
%     vo.y_names = {'share'};
%     vo.state_names = {'Qshare', 'PE'};
% elseif strcmpi(vo.model,'regret')
%     vo.evo_fname = @f_trust_Qlearn_counter_hybrid_regret_pmv;
%     vo.hidden_states=2;
%     vo.n_theta=1; %learning rate
%     vo.theta_names={'alpha'};
%     vo.state_names = {'Qshare', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model,'null')
%     vo.evo_fname = @f_trust_Qlearn_null_pmv;
%     vo.hidden_states=2;
%     vo.n_theta=1; %learning rate alpha and uncertainty sensitivity parameter tau
%     vo.theta_names={'alpha'};
%     vo.state_names = {'Qshare', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'countertrustee')
%     vo.evo_fname = @f_trust_Qlearn_counter_trustee;
%     vo.hidden_states=2; %track basis weights (value), as well as PE as tag-along state
%     vo.n_theta=1; %learning rate and selective maintenance
%     vo.theta_names={'alpha'};
%     vo.state_names = {'Qshare', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model,'ptrustee')
%     vo.evo_fname = @f_trust_Qlearn_ptrustee;
%     vo.hidden_states = 2;
%     vo.n_theta = 3;
%     vo.theta_names = {'alpha_bad', 'alpha_neutral', 'alpha_good'};
%     vo.state_names = {'Qshare', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model,'prs')
%     vo.evo_fname = @f_trust_Qlearn_prs;
%     vo.hidden_states = 2;
%     vo.n_theta = 3;
%     vo.theta_names = {'alpha_poor', 'alpha_neutral', 'alpha_rich'};
%     vo.state_names = {'Qshare', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model,'discrepancy')
%     vo.evo_fname = @f_trust_Qlearn_subjective_discrepancy;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.n_theta = 5;
%     vo.theta_names = {'alpha', 'cooperate', 'defect', 'relief', 'regret'};
%     vo.n_outputs = 1; %two actions + no response
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'twoQs_Actual_Outcome')
%     vo.evo_fname = @f_trust_Qlearn_twoQs_Actual_Outcome;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'twoQs_Counterfactual')
%     vo.evo_fname = @f_trust_Qlearn_twoQs_Counterfactual;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'twoQs_Countertrustee')
%     vo.evo_fname = @f_trust_Qlearn_twoQs_Countertrustee;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'twoQs_Full_Regret')
%     vo.evo_fname = @f_trust_Qlearn_twoQs_Full_Regret;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'twoQs_Original_Regret')
%     vo.evo_fname = @f_trust_Qlearn_twoQs_Original_Regret;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.y_names = {'share'};
% elseif strcmpi(vo.model, 'twoQs_Policy')
%     vo.evo_fname = @f_trust_Qlearn_twoQs_Policy;
%     vo.hidden_states = 3;
%     vo.state_names = {'Qshare', 'Qkeep', 'PE'};
%     vo.y_names = {'share'};
% end

% 
%     if vo.sigmakappa == 0
%         vo.obs_fname=@g_trust_softmax1;
%         vo.phi_names = {'beta'};% observation function (softmax mapping), evaluates Q(share)
%     elseif vo.sigmakappa == 1 && counter == 0
%         vo.obs_fname = @g_trust_softmax_ED1;% observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter, w/ choice rule appropriate for actual evolution function
%         vo.phi_names = {'beta', 'kappaS'};
%     elseif vo.sigmakappa == 1
%         vo.obs_fname = @g_trust_softmax_ED2;  % observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter, w/ choice rule appropriate for other than actual evolution function(s)
%         vo.phi_names = {'beta', 'kappaS'};
%     elseif vo.sigmakappa == 2
%         vo.obs_fname = @g_trust_softmax_ks_kt;% observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter + trustee-wise kappa parameter (requires multisession)
%         vo.phi_names = {'beta', 'kappaS', 'kappaST'};
%     elseif vo.sigmakappa ==3
%         vo.obs_fname = @g_trust_softmax_ks_kt_factorized;
%         vo.phi_names = {'beta', 'kappaS', 'kappaST'};
%     elseif vo.sigmakappa ==4
%         vo.obs_fname = @g_trust_softmax_ks_kt_nomult_dummy;
%         vo.phi_names = {'beta', 'kappaS', 'kappaSTD1', 'kappaSTD2'};
% 
%     elseif vo.sigmakappa == 5
%         if strcmpi(vo.dataset,'learn') || vo.fourTrustee
%             vo.obs_fname = @g_trust_softmax_ks_kt_nomult_dummy_fourT;
%             vo.phi_names = {'beta', 'kappaS', 'kappaSTD1', 'kappaSTD2','kappaSTD3'};
%         else
%             vo.obs_fname = @g_trust_softmax_ks_kt_nomult_dummy_fixed0;
%             vo.phi_names = {'beta', 'kappaS'};
%         end
%     end
% 
% 
% if vo.sigmakappa > 0
%     if vo.sigmakappa < 3
%         vo.n_phi = 1+vo.sigmakappa;%observation function parameters: by default = 1 (beta), +1 kappa for the subject's altruistic bias, +1 for the subject's trustee specific bias
%     elseif vo.sigmakappa == 5 && ~(strcmpi(vo.dataset,'learn') || vo.fourTrustee)
%         vo.n_phi = 2;
%     else
%         vo.n_phi = vo.sigmakappa;
%     end
% else
%     vo.n_phi = 1;
% end
% 
% 
% if strcmpi(vo.model,'discrepancy')
%     vo.n_phi = 4;
%     vo.obs_fname=@g_trust_two_Qs;
%     vo.phi_names = {'beta', 'kappaS', 'kappaSTD1', 'kappaSTD2'};
% elseif strcmpi(vo.model, 'rationalsd')
%     vo.obs_fname = @g_trust_softmax_ks_kt_fivePhi;
%     vo.n_theta=1;
%     vo.n_phi=5;
%     vo.phi_names={'beta', 'kappaS', 'kappaST', 'share_beta', 'keep_beta'};
%     vo.theta_names={'alpha'};
% % elseif strcmpi(vo.model, 'twoQs_Actual_Outcome') | strcmpi(vo.model, 'twoQs_Counterfactual') | strcmpi(vo.model, 'twoQs_Countertrustee') | strcmpi(vo.model, 'twoQs_Full_Regret') | strcmpi(vo.model, 'twoQs_Original_Regret') | strcmpi(vo.model, 'twoQs_Policy')
% %     vo.obs_fname = @g_trust_two_Qs;
% %     vo.n_theta=1;
% %     vo.n_phi=3;
% %     vo.phi_names={'beta', 'kappaS', 'kappaST'};
% %     vo.theta_names={'alpha'};
% elseif strcmpi(vo.model, 'ratrightness') | strcmpi(vo.model, 'rationalpol') | strcmpi(vo.model,'rationalreg') | strcmpi(vo.model,'rationalremod') | strcmpi(vo.model,'rationalcounter') | strcmpi(vo.model,'rationalfullpol') | strcmpi(vo.model,'utilitygain') 
%     vo.obs_fname = @g_trust_softmax_ks_kt;
%     vo.n_theta=1;
%     vo.n_phi=3;
%     vo.phi_names={'beta', 'kappaS', 'kappaST'};
%     vo.theta_names={'alpha'};
% elseif strcmpi(vo.model, 'ratcoop')
%     vo.obs_fname = @g_trust_softmax_ks_kt_four_phi;
%     vo.n_theta=1;
%     vo.n_phi=4;
%     vo.phi_names={'beta', 'kappaS', 'kappaST', 'omega'};
%     vo.theta_names={'alpha'};
% end

end
