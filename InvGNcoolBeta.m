function [mpred, r, info] = InvGNcoolBeta(funFwd, dobs, m0, target, options, params)

% initialize outputs
r     = nan;
mpred = m0;
info  = struct;

%% Check for parameters, if not there, set defaults

% max number of iterations
if ~isfield(options,'iterations')
    maxIt = 100;
else
    maxIt = options.iterations;
end

if ~isfield(options,'iterationsLS')
    maxItLS = 20;
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
    LSshorten = 0.25;
else
    LSshorten = options.LSshorten;
end

if ~isfield(options,'tolF')
    tolF = 1e-1;
else
    tolF = options.tolF;
end

if ~isfield(options,'tolX')
    tolX = 1e-16;
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

%% Set Beta schedule
% if no initial beta given, estimate it by one iteration of the power
% method (with magic number 100!)

if ~isfield(params,'beta0')
    x0 = rand(numel(m0),1);
    [~,J0] = funFwd(m0);
    top = norm(Wd*(J0*x0),2);
    bot = norm(Wm*x0,2);
    beta0 = 10^2*sqrt(top/bot);
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

%%
% initialize space for data misfit, regularization and objective function
fphid = zeros(maxIt,1);
fphim = zeros(maxIt,1);
f     = zeros(maxIt,1);
beta  = zeros(maxIt,1);
normg = zeros(maxIt,1);

% create schedule
betaupdate          = false(maxIt,1);
betaind             = 1:coolingrate:maxIt;
betaupdate(betaind) = true;

% start log
logD = ' %5i  %8i %9.4e  %9.4e  %9.2e  %9.1f %9.1f\n';
logH = ' %5s  %8s %9s  %9s  %9s  %9s %9s \n';
fprintf('\n\n -------------- INVERSION W/ GN & BETA COOLING ------------- \n');
fprintf(logH,'Iter','LS Iter','PhiD','PhiM','Beta','Objective','Norm g');

% initialize
m        = m0;
[dpred,~]   = funFwd(m);
r        = dobs - dpred;
fphid(1) = funPenalty(r,params);
fphim(1) = primal_norm(m,Wm,params);

beta(1) = beta0;
f(1) = fphid(1) + beta(1)*fphim(1);

LSit = 0; 
%print
%fprintf(logD,1,LSit,fphid(1),fphim(1),beta(1),phi(1),normg(1));

%% perform inversion
for i = 2:maxIt
    
    % test if we have reached target misfit
    if fphid(i-1) <= target
        mpred = m;
        info.fphid = fphid(1:i);
        info.fphim = fphim(1:i);
        info.phi   = f(1:i);
        info.beta  = beta(1:i);
        info.normg = normg(1:i);
        info.exit  = 'Target';
        fprintf(' --------- Target Reached ---------- \n');
        return
    end
    
    % check if beta should be update
    if betaupdate(i)
        beta(i) = beta(i-1)*coolingfact;
    else
        beta(i) = beta(i-1); 
    end
    
    % compute data, resudial
    [dpred,J] = funFwd(m);
    r         = dobs-dpred;
    dr        = -1; 
    
    % compute phid, phim
    [fphid(i),gphid,Hphid] = funPenalty(r,params);
    [fphim(i),gphim,Hphim] = primal_norm(m,Wm,params);
    
    f(i) = fphid(i) + beta(i)*fphim(i);
    g    = dr*J'*gphid + beta(i)*gphim;
    H    = J'*(Hphid)*J + beta(i)*Hphim; %J'*(Hphid + gphid*gphid')*J + beta(i)*Hphim;

    normg(i) = norm(g);
    
    % print log
    fprintf(logD,i-1,LSit,fphid(i),fphim(i),beta(i),f(i),normg(i));
    
    if normg(i) < tolG
        fprintf('\n ---- |g| < tolG ------ \n');
        mpred = m;
        info.exit = 'tolG';
        return
    end
    
    % compute step
    step = -H\g;
    clear H

    % line search
    LSit = 0;

    while 1
        if LSit > maxItLS % cool beta and try again
            mpred = m; 
            beta(i) = beta(i-1)*coolingfact;
            betaupdate(i:end) = false;
            betaupdate(i:coolingrate:end) = true;
            warning('LineSearch:maxItLS','Line search broke :(');
            info.exit = 'maxItLS';
            break
        end
        
        % see if there is sufficient change in the model - if not, cool
        % beta and start again
        if norm(step) < tolX
            warning('Linesearch:normX','Insufficient change in model, beta reduced'); 
            beta(i) = beta(i-1)*coolingfact;
            betaupdate(i:end) = false;
            betaupdate(i:coolingrate:end) = true;
            info.exit = 'tolX';
            break
        end
        
        % check if step is too large, if so... scale back
        if max(abs(step)) > maxStep
            step = step * maxStep/max(abs(step));
            fprintf('   ---- Step Scaled back ---- \n'); 
        end
        
        % take a step
        mtry = m+step;
 
        
        % compute data, resudial
        dtry = funFwd(mtry);
        rtry = dobs-dtry;
        
        % compute phid, phim, ftry
        [fphidtry] = funPenalty(rtry,params);
        [fphimtry] = primal_norm(mtry,Wm,params);
        ftry       = fphidtry + beta(i)*fphimtry;
        
        % LSshorten
        if f(i) - ftry > tolF % if sufficient decrease, exit line search
            break
        end
        step = step*LSshorten;
        LSit = LSit + 1;
    end
    fphid(i) = fphidtry;
    fphim(i) = fphimtry; 
    f(i)     = ftry; 
    m = mtry; % update model
    
    
    
end

mpred = m;
warning('GN:maxit','Maximum number of GN iterations reached');
info.exit = 'maxIT';


end