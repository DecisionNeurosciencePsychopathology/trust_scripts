function [options, dim] = mmrescan_get_vba_options(so)
%% this function sets up the options and dim structures for VBA

so=mmrescan_validate_options(so); %should be handled upstream, but just in case

options=[];
priors=[];

if ~so.graphics
  options.DisplayWin = 0; %whether to display graphics during fitting
  options.GnFigs = 0;
end

% u is 2 x ntrials where first row is rt and second row is reward

% copy so to inF and inG as starting point
options.inF = so;
options.inG = so;
% convergence settings
options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=0; %don't show single subject fitting process
%% split into conditions/runs
if so.multisession
    options.multisession.split = repmat(so.n_t/so.n_trustees,1,so.n_trustees); %splitting the sessions
    %% fix parameters
    if so.fixed == 1
        %options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        %options.multisession.fixed.phi = 'all';
        options.multisession.fixed.theta = 'all';
        %options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions
        options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        %options.multisession.fixed.X0 = 'all';
        %% set priors on initial states
        priors.SigmaX0 = diag([.3 0]);  %X0 is allowed to vary between sessions
    elseif so.fixed == 2
        %options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        %options.multisession.fixed.phi = 'all';
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions; subject-wise kappa varies between sessions
        %options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        options.multisession.fixed.X0 = 'all';
        priors.SigmaX0 = diag([0 0]);   %infinite precision priors set on the initial value and PEs
    elseif so.fixed == 3
        %options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        %options.multisession.fixed.phi = 'all';
        options.multisession.fixed.theta = 'all';
        %options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions
        options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        options.multisession.fixed.X0 = 'all';     
        priors.SigmaX0 = diag([0 0]);   %infinite precision priors set on the initial value and PEs
    elseif so.fixed == 4
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        options.multisession.fixed.X0 = 'all';
        priors.muX0 = diag([0 0]);
        priors.SigmaX0 = diag([0 0]);
    
    end        
end
skipf = zeros(1,so.n_t);
skipf(1) = 1; %skip first trial
options.skipf = skipf;
options.isYout = zeros(1,so.n_t);
%% specify dimensions of data to be fit
dim = struct('n', so.hidden_states, ...
  'n_theta', so.n_theta, ...
  'n_phi', so.n_phi, ...
  'n_t', so.n_t, ...
  'p', 1); %isYout should be 1* n_t

%%populate priors
priors = mmrescan_get_priors(dim, so);

options.priors = priors;
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)

end
