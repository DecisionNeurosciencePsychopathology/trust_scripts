function [so] = mmrescan_validate_options(so)
% nbasis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
% fixed_params_across_runs -- self-explanatory
% fit_propspread -- makes temporal generalization within the eligibility trace a free parameter
% ntimesteps:      number of time bins
% u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum

if nargin < 1, so=[]; end

%if user specifies a dataset upstream, don't check environment variable
if isfield(so, 'dataset')
    fprintf('Using user-specified dataset: %s.\n   Will not check dataset environment variable', so.dataset);
else
    so.dataset=getenv('dataset');
    if strcmpi(so.dataset, '')
        so.dataset='mmrescan';
    end
end

%same principle for model variant
if isfield(so, 'model')
    fprintf('\nUsing user-specified model: %s.\n   Will not check model environment variable', so.model);
else
    so.model=getenv('model');
    if strcmpi(so.model, '')
        so.model = 'policy'; % will run to get value and prediction errors.
    end
end

if ~isfield(so, 'multisession'), so.multisession=0; end
if ~isfield(so, 'graphics'), so.graphics=0; end
if ~isfield(so, 'model'), so.model='policy'; end
if ~isfield(so, 'n_trials'), so.n_trials=144; end
if ~isfield(so, 'saveresults'), so.saveresults=0; end %used in some places to denote whether to save fitted outputs
if ~isfield(so, 'params'), so.params = 'notfixedparams'; end

%hidden states field used to index hidden state vector inside evolution and observation functions
if strcmpi(so.model,'policy')
    so.evo_fname = @f_trust_Qlearn_policy;
    so.hidden_states=2; %track basis weights (value), as well as PE as tag-along state
    so.n_theta=1; %learning rate
    so.theta_names={'alpha'};
elseif strcmpi(so.model, 'mixed')
    so.evo_fname = @f_trust_Qlearn_policy_actual_mix;
    so.hidden_states = 2;
    so.n_theta=2; %learning rate and mixing parameter for the compound value
    so.theta_names={'alpha','mix_param'};
elseif strcmpi(so.model, 'mixed2')
    so.evo_fname = @f_trust_Qlearn_policy_actual_mix2;
    so.hidden_states = 2;
    so.n_theta=2; %learning rate and mixing parameters for the compound value
    so.theta_names={'alpha','mix_param'};
elseif strcmpi(so.model, 'rcounter') || strcmp(so.model, 'rcounter2') || strcmp(so.model,'rcounter3')
    so.evo_fname = @f_trust_Qlearn_rcounter;
    so.hidden_states=2; %track basis weights (value), as well as PE as tag-along state
    so.n_theta=1; %learning rate
    so.theta_names={'alpha'};
elseif strcmpi(so.model,'utilitygain')
    so.evo_fname = @f_trust_Qlearn_utilitygain;
    so.hidden_states=2;
    if so.utility_goodbad
        so.n_theta=2;
        so.theta_names={'alpha_good','alpha_bad'};
    else
        so.n_theta=1;
        so.theta_names={'alpha'};
    end
elseif strcmpi(so.model,'regret')
    so.evo_fname = @f_trust_Qlearn_counter_hybrid_regret_pmv;
    so.hidden_states=2;
    so.n_theta=1; %learning rate
    so.theta_names={'alpha'};
elseif strcmpi(so.model,'null')
    so.evo_fname = @f_trust_Qlearn_null_pmv;
    so.hidden_states=2;
    so.n_theta=1; %learning rate alpha and uncertainty sensitivity parameter tau
    so.theta_names={'alpha'};
elseif strcmpi(so.model, 'countertrustee')
    so.evo_fname = @f_trust_Qlearn_counter_trustee;
    so.hidden_states=2; %track basis weights (value), as well as PE as tag-along state
    so.n_theta=1; %learning rate and selective maintenance
    so.theta_names={'alpha'};
elseif strcmpi(so.model,'ptrustee')
    so.evo_fname = @f_trust_Qlearn_ptrustee;
    so.hidden_states = 2;
    so.n_theta = 3;
    so.theta_names = {'alpha_bad', 'alpha_neutral', 'alpha_good'};
elseif strcmpi(so.model,'prs')
    so.evo_fname = @f_trust_Qlearn_prs;
    so.hidden_states = 2;
    so.n_theta = 3;
    so.theta_names = {'alpha_poor', 'alpha_neutral', 'alpha_rich'};
elseif strcmpi(so.model,'discrepancy')
    so.evo_fname = @f_trust_Qlearn_subjective_discrepancy;
    so.hidden_states = 3;
    so.n_theta = 5;
    so.theta_names = {'alpha', 'cooperate', 'defect', 'relief', 'regret'};
end


if strcmpi(so.model,'discrepancy')
    so.obs_fname=@g_trust_two_Qs;
    so.phi_names = {'beta'};
else
    if so.sigmakappa == 0
        so.obs_fname=@g_trust_softmax1;
        so.phi_names = {'beta'};% observation function (softmax mapping), evaluates Q(share)
    elseif so.sigmakappa == 1 && counter == 0
        so.obs_fname = @g_trust_softmax_ED1;% observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter, w/ choice rule appropriate for actual evolution function
        so.phi_names = {'beta', 'kappaS'};
    elseif so.sigmakappa == 1
        so.obs_fname = @g_trust_softmax_ED2;  % observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter, w/ choice rule appropriate for other than actual evolution function(s)
        so.phi_names = {'beta', 'kappaS'};
    elseif so.sigmakappa == 2
        so.obs_fname = @g_trust_softmax_ks_kt;% observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter + trustee-wise kappa parameter (requires multisession)
        so.phi_names = {'beta', 'kappaS', 'kappaST'};
    elseif so.sigmakappa ==3
        so.obs_fname = @g_trust_softmax_ks_kt_factorized;
        so.phi_names = {'beta', 'kappaS', 'kappaST'};
    elseif so.sigmakappa ==4
        so.obs_fname = @g_trust_softmax_ks_kt_nomult_dummy;
        so.phi_names = {'beta', 'kappaS', 'kappaSTD1', 'kappaSTD2'};
    elseif so.sigmakappa == 5
        if strcmpi(so.dataset,'learn') || so.fourTrustee
            so.obs_fname = @g_trust_softmax_ks_kt_nomult_dummy_fourT;
            so.phi_names = {'beta', 'kappaS', 'kappaSTD1', 'kappaSTD2','kappaSTD3'};
        else
            so.obs_fname = @g_trust_softmax_ks_kt_nomult_dummy_fixed0;
            so.phi_names = {'beta', 'kappaS'};
        end
    end
end
if so.sigmakappa > 0
    if so.sigmakappa < 3
        so.n_phi = 1+so.sigmakappa;%observation function parameters: by default = 1 (beta), +1 kappa for the subject's altruistic bias, +1 for the subject's trustee specific bias
    elseif so.sigmakappa == 5 && ~(strcmpi(so.dataset,'learn') || so.fourTrustee)
        so.n_phi = 2;
    else
        so.n_phi = so.sigmakappa;
    end
else
    so.n_phi = 1;
end
if strcmpi(so.model,'discrepancy')
    so.n_phi = 4;
    so.phi_names = {'beta', 'kappaS', 'kappaSTD1', 'kappaSTD2'};
end

so.trials_per_run = 48;
%so.n_t = so.n_trials;
so.n_hidden_states = so.hidden_states;
end
