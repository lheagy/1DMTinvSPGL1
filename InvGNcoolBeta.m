function [mpred, r, info] = InvGNcoolBeta(funFwd, dobs, target, m0, options, params)

%% Check for parameters, if not there, set defaults

% max number of iterations
if ~isfield(options,'iterations')
    maxIt = 100;
else
    maxIt = options.iterations;
end

if ~isfield(options,'iterationsLS')
    maxItLS = 10;
else
    maxItLS = options.iterationsLS;
end

if ~isfield(options,'maxStep')
    maxStep = 10^6;
else
    maxStep = options.maxStep;
end    

if ~isfield(options,'LSreduction');
    LSreduction = 1e-4;
else 
    LSreduction = options.LSreduction;
end

if ~isfield(options,'LSshorten')
    LSshorten = 0.5;
else
    LSshorten = options.LSshorten;
end

if ~isfield(options,'tolF')
    tolF = 1e-1;
else 
    tolF = options.tolF;
end

if ~isfield(options,'tolX')
    tolX = 1e-1;
else
    tolX = options.tolX;
end

if ~isfield(options,'tolG')
    tolG = 1e-1;
else
    tolG = options.tolG;
end

%% Grab model and data weights
if ~isfield(params,'Wd')
    Wd = speye(numel(f));
else
    Wd = params.Wd;
end

if ~isfield(params,'Wm')
    Wm = speye(numel(m0));
else
    Wm = params.Wm;
end

% function handles for data misfit and model regularization
funPenalty  = options.funPenalty;
primal_norm = options.primal_norm; 
primal_normDeriv = options.primal_normDeriv;

%% Set Beta schedule
% if no initial beta given, estimate it by one iteration of the power
% method (with magic number 100!) 
if ~isfield(params,'beta0')
    x0 = rand(numel(m0),1);
    [~,J0] = funFwd(m0);
    top = x0'*(J0'*(Wd'*(Wd*(J0*x0))));
    bot = x0'*(Wm'*(Wm*x0)); 
    beta0 = 100*sqrt(top/bot);
else
    beta0 = params.beta0;
end

% beta cooling schedule
if ~isfield(params,'coolingrate')
    coolingrate = 2;
else
    coolingrate = params.coolingrate;
end

% beta cooling factor
if ~isfield(params,'coolingfact');
    coolingfact = 1/2;
else
    coolingfact = params.coolingfact;
end

% initialize space for data misfit, regularization and objective function
phid = zeros(maxIt+1,1);
phim = zeros(maxIt+1,1);
phi  = zeros(maxIt+1,1);
beta = zeros(maxIt+1,1); 

% create schedule
betaupdate = false(maxIt,1);
betaind = 1:coolingrate:maxIt;
betaupdate(betaind) = true;

% start log
logD = ' %5i  %9.4e  %9.4e  %9.2e  %9.1f \n';
logH = ' %5s  %10s  %10s  %9s  %9s \n';
fprintf('\n\n --------- INVERSION W/ GN & BETA COOLING --------- \n'); 
fprintf(logH,'Iter','PhiD','PhiM','Beta','Objective');

% initialize
[d0,~] = funFwd(m0);
r0      = dobs - d0;
phid(1) = funPenalty(r0,params); 

phim(1) = primal_norm(m0,Wm,params);

beta(1) = beta0;
phi(1) = phid(1) + beta(1)*phim(1);

%print
fprintf(logD,0,phid(1),phim(1),beta(1),phi(1));

%% perform inversion
for i = 1:maxIt
    % test if we have reached target misfit
    if phid(i+1) <= target
        mpred = m0;
        break
    end
    
    dpred = funFwd(m0);
    r     = dobs-dpred;
    [phid(i+1),dphid] = phid(r,params);
    []
    
    % check if beta should be update
    if betaupdate(i) && i > 1
        beta(i) = beta(i-1)*betafact;
    end
    
    % compute GN update
    
    
    
    fprintf(logD,0,phid(i+1),phim(i+1),beta(i+1),phi(i+1));
    
    
end



end