function [priors] = mmrescan_get_priors(dim, so)
priors=[];
sigmakappa = so.sigmakappa;
multisession = so.multisession;
fixed = so.fixed;
if dim.n_phi == 2
    priors.muPhi = [1;0];
elseif dim.n_phi == 3
    priors.muPhi = [1;0;0];
elseif dim.n_phi == 4
    priors.muPhi = [1;0;0;0];
else
    priors.muPhi = 1;
end

    
%setting 0 mean prior on Q values
priors.muX0 = zeros(so.n_hidden_states,1);

if strcmpi(so.model,'discrepancy')
    priors.muTheta = zeros(dim.n_theta,1); 
    priors.SigmaTheta = diag(ones(1,dim.n_theta));
else
    priors.muTheta = zeros(dim.n_theta,1); 
    priors.SigmaTheta = diag(10*ones(1,dim.n_theta));
end


    


if sigmakappa > 0
    priors.SigmaPhi = diag([10 1/3*ones(1,dim.n_phi-1)]); %better than when xcompared with sigma_kappa = 1
    %priors.SigmaPhi = diag([10 ones(1,dim.n_phi-1)]);
else
    priors.SigmaPhi = 10;
end


priors.a_alpha = Inf;   % infinite precision prior
priors.b_alpha = 0;
priors.a_sigma = 1;     % Jeffrey's prior
priors.b_sigma = 1;     % Jeffrey's prior

%0 mean and variance on initial states
priors.muX0 = zeros(dim.n,1);
priors.SigmaX0 = zeros(dim.n);

if multisession
    options.multisession.split = repmat(so.n_t/so.n_trustees,1,so.n_trustees); %splitting the sessions
    % fix parameters
    if fixed == 1
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        priors.SigmaX0 = diag([.3 0]);  %X0 is allowed to vary between sessions
    elseif fixed == 2
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions; subject-wise kappa varies between sessions
        options.multisession.fixed.X0 = 'all';
        priors.SigmaX0 = diag([0 0]);   %infinite precision priors set on the initial value and PEs
    elseif fixed == 3
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        options.multisession.fixed.X0 = 'all';
        priors.SigmaX0 = diag([0 0]);   %infinite precision priors set on the initial value and PEs
    end
end

% if strcmpi(so.params, 'fixedparams')
%     dim.n_phi =4;
%     dim.n_theta = 1;
%     priors.muPhi = [-0.2450279; -0.4439473; -0.2729371;-0.006007666];
%     priors.SigmaPhi = zeros(dim.n_phi);
%     priors.muTheta = [-0.3562523];
%     priors.SigmaTheta = zeros(dim.n_theta);
% end

