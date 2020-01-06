function  [ gx] = g_trust_softmax_ks_kt_2actions(x,phi,u,in )
% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = exp(phi(1));
kappaS = phi(2);  %subject-wise bias parameter
kappaST = phi(3); %trustee-wise bias parameter

dQ = (x(1)-x(2));

gx = sig(kappaS + kappaST + beta*dQ); %edited observation function with a subject-specific bias parameter + subject specific trustee-specific bias parameter

dgdx = zeros(size(x,1),1);
dgdx(1) = beta*gx*(1-gx);
dgdx(2) = -beta*gx*(1-gx);
dgdphi = [beta*dQ*gx*(1-gx),gx*(1-gx)];

% dgdx = zeros(size(x,1),1); % for value tracking only
% dgdx = zeros(size(1,1),1);  % for value + pe tracking
% dgdx(1) = beta*gx*(1-gx);
% dgdP(1) = beta*dQ*gx*(1-gx);
