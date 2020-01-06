function  [fx] = f_trust_Qlearn_policy_actual_mix(x,theta,u,inF)

% function  [fx,dfdx,dfdP] = f_trust_Qlearn_counter(x,theta,u,inF)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : q-values (2x1)
%   - P : (inverse-sigmoid) learning-rate
%   - u : u(1)=previous action (1 or 0), u(2)=feedback
%   - in : [useless]
% OUT:
%   - fx: evolved q-values (2x1)
%   - dfdx/dfdP: gradient of the q-values evolution function, wrt q-values
%   and evolution parameter, respectively.F 

% u(1) = subject's decisions
% u(2) = trustee's decisions


alpha = 1./(1+exp(-theta(1))); % learning rate is bounded between 0 and 1.
p1 = 1./(1+exp(-theta(2))); %parameter regulating contribution of actual rewards to the total, also bounded between 0 and 1.
p2 = 1 - p1; %parameter regulating contribution of policy-style rewards to the total

%reward for trial 
if (u(2)==1 && u(1)==1)     %trustee shared, subject shared
    r_a = 1.5;              
    r_p = 1;
elseif (u(2)<1 && u(1)==1) %trustee kept, subject shared
    r_a = 0;
    r_p = 0;
elseif (u(2)<1 && u(1)<1)  %trustee kept, subject kept
    r_a = -1; %was 1
    r_p = -1; %was 1
else
    r_a = -1; %was 1              %trustee shared, subject kept
    r_p = 0;
end

r = p1 * r_a + p2 * r_p;

fx = zeros(length(x),1);
pe = r-x(1); % prediction error for share
fx(1) = x(1) + alpha*pe;
fx(2) = pe;

%% one hidden state (value)
% dfdx = zeros(size(x,1),1);
% dfdx(1) = [1-alpha];
% dfdP = [alpha*(1-alpha)*pe];

%% two hidden states (value + pe)
% gradients' derivation
% if u(1)==1
%     dfdx = [1-alpha, 0;
%             0, 1];
%     dfdP = [alpha*(1-alpha)*pe(1),0];
% else
%     dfdx = [1, 0;
%             0, 1-alpha];
%     dfdP = [0,alpha*(1-alpha)*pe(2)];
% end
