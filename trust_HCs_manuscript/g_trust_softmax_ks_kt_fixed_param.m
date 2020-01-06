function  [ gx] = g_trust_softmax_ks_kt_fixed_param(x,phi,u,in )
% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : 
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = exp(phi(1));
kappaS = phi(2);  %subject-wise bias parameter
%kappaST = phi(3); %trustee-wise bias parameter

if u(8) == 1 %good trustee
    kappaST = 0.03307049;
elseif u(8) == 2 %bad trustee
    kappaST = 0.0102719;
elseif u(8) == 3 %neutral trustee
    kappaST = -.00496066;
else %computer trustee
    kappaST = -0.02498429;
end

dQ = x(1);
%dQ = x(1)-1; %for the corrected null model

gx = sig(kappaS + kappaST + beta*dQ); %edited observation function with a subject-specific bias parameter + subject specific trustee-specific bias parameter

% dgdx = zeros(size(x,1),1); % for value tracking only
% dgdx = zeros(size(1,1),1);  % for value + pe tracking
% dgdx(1) = beta*gx*(1-gx);
% dgdP(1) = beta*dQ*gx*(1-gx);
