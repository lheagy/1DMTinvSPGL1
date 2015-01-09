function [f,g,H] = objectiveTradeOff(dobs,m,phidfun,phimfun,funFwd,params)

if ~isfield(params,'beta')
    warning('Beta:notDefined','Beta not defined, assigned value of 1');
    beta = 1;
else
    beta = params.beta;
end

Wm   = params.Wm; 

% compute data, resudial
[dpred,J] = funFwd(m);
r     = dobs-dpred;

% compute phid, phim
[fphid,gphid,Hphid] = phidfun(r,params);
[fphim,gphim,Hphim] = phimfun(m,Wm,params);

% compute f, g, for full objective function
f = fphid + beta*fphim;
g = J'*gphid + beta*gphim;
H = J'*(Hphid + gphid'*gphid)*J + beta*Hphim;


end