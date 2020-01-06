function [theta_trans] = transform_theta(theta, inF)

if any(strcmpi(inF.model, {'fixed', 'decay', 'fixed_UV','utilitygain', 'decay_factorize', 'policy', 'countertrustee', 'null', 'regret', 'ptrustee', 'prs','discrepancy','mixed','rcounter','rcounter2','rcounter3', 'mixed2'}))
  theta_trans = VBA_sigmoid(theta); %sigmoid transform alpha, which is theta(1)
else
  error(['unrecognized model in transform_theta: ', inF.model]);
end

%handle other models at some point...

end
