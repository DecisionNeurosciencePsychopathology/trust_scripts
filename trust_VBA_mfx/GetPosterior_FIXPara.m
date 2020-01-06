function [posterior,out]=GetPosterior_FIXPara(u,y,options,muPhi,muTheta)

    dim = options.dim;
    
    options.priors.muPhi = muPhi;
    options.priors.muTheta = muTheta;
    
    options.priors.SigmaPhi = zeros(length(muPhi),length(muPhi));
    options.priors.SigmaTheta = zeros(length(muTheta),length(muTheta));
    
    options.priors.SigmaX0 = zeros(dim.n);
    
    options.DisplayWin = 0;
    %options.isYout = zeros(length(y));
    [posterior,out] = VBA_NLStateSpaceModel(y,u,options.f_fname,options.g_fname,dim,options);

end