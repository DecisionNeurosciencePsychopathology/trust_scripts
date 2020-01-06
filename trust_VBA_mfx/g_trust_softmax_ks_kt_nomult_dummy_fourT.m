function  [ gx] = g_trust_softmax_ks_kt_nomult_dummyLearn(x,phi,u,inG )
% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = exp(phi(1));
kappaS = phi(2);  %subject-wise bias parameter
%kappaST = phi(3); %trustee-wise bias parameter

trustee=u(3);
if trustee == 0 %neutral
    kappaST = 0; 
elseif trustee == 1 %computer
    kappaST = phi(3); 
elseif trustee == 2 %good
    kappaST = phi(4);
elseif trustee == 3 %bad
    kappaST = phi(5);
else
    error('Unknown trustee code');
end
    
dQ = x(1);

if strcmpi(inG.model,'utilitygain') || strcmpi(inG.model,'rcounter2')
    dQ = (1.5 * dQ) - 1;
elseif strcmp(inG.model, 'rcounter')
    dQ = 1.5 * dQ;
end
%dQ = x(1)-1; %for the corrected null model

gx = sigmoid(kappaS + kappaST + beta*dQ); %edited observation function with a subject-specific bias parameter + subject specific trustee-specific bias parameter

% dgdx = zeros(size(x,1),1); % for value tracking only
% dgdx = zeros(size(1,1),1);  % for value + pe tracking
% dgdx(1) = beta*gx*(1-gx);
% dgdP(1) = beta*dQ*gx*(1-gx);
