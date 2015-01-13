% 1DFDEM Inversion
% 
% Driver file for the 1D magnetotelluric forward and inverse problems.
%
%
% Lindsey J. Heagy
% lheagy@eos.ubc.ca
%
% ------------------------------------------------------------------------%

%% House-keeping 
clear all
close all
% clc

testit = false; 

% add path to SPGL1general
addpath(genpath('../spgl1'));
if testit; addpath('../tests/'); ntest = 7; end

%% Parameters

% mesh
nc   = 91;      % number of cells

% frequency 
frange = [-2, 2]; % log10 of fmin, fmax
nf   = 10;      % number of frequencies 
nd   = nf*2;    % number of data

% reference model
mref = -1;      % log of reference sigma
mu0  = 4e-7*pi; % magnetic permeability

%% Mesh
mesh = getMesh(nc,exp(mref),10.^frange,exp(-1),3e5);

%% Model

eta = log(10^4);
%t = log(mesh.zc);
t = log10(mesh.zc); 
mb = -1;
%mtrue = mb*ones(mesh.nzc,1);
%mtrue((t > 0.5)&(t < 2)) = 1;
%mtrue((t > 2.5)&(t< 3.5)) = 0.75; 
mtrue = mb+1*exp(-(t-0.25*eta).^2); %+1.5*exp(-(t-0.85*eta).^2/4);
mtrue = mtrue(:); 
model.mu = mu0*ones(mesh.nzc,1);

model.transform = @(m) exp(m);
model.itransform = @(sigma) log(sigma);
model.transformderiv = @(m) sdiag(exp(m));
model.mref = mref;

m = mtrue-mref;
model.m = m;

% if testit
%     fprintf('Testing Model Mapping \n'); 
%     fun = @(m) {model.transform(m),model.transformderiv(m)};
%     testGrad(fun,exp(rand(numel(m),1)));
%     fprintf('\n\n');
% end
%% forward model

f = logspace(frange(1),frange(2),nf);
dtrue = get1DMTfwd(m,model,mesh,f);

%% Set up weights for objective functions

nzc = mesh.nzc; 

% add noise to data
percn = 0.03;
noise = percn*randn(nd,1).*abs(dtrue);
dobs = dtrue+noise;

% Data weight
perc = 0.05;
flr  = 0.003;
flrr  = flr*mean(abs(dobs(1:nf)));
flri  = flr*mean(abs(dobs(nf+1:end)));
Wd   = perc*abs(dobs) + [flrr*ones(nf,1);flri*ones(nf,1)];
Wd   = sdiag(1./Wd);

% Model Weights
alphas  = 1.5;
alphaz  = 3;
alpha2z = 10; 
Wms = speye(mesh.nzc);

t   = log(mesh.z);
dt  = log(mesh.dz);
dts = (dt(1:end-1) + dt(2:end))/2;
Wmd = sdiag([1./dts])*spdiags([-1*ones(nzc,1), ones(nzc,1)],-1:0,nzc,nzc);
Wmd(1,1) = 0;

% 
Wmd2 = sdiag([1./dts]).^2*spdiags([-1*ones(nzc,1), 2*ones(nzc,1), -1*ones(nzc,1)],[-1:1],nzc,nzc);
Wmd2(1,[1:2]) = 0; 
Wmd2(end,[end-1:end]) = 0; 

Wm  = [alphas*Wms; alphaz*Wmd; alpha2z*Wmd2];
[nWm,mWm] = size(Wm);

if nWm~=mWm
    Wm2 = Wm'*Wm;
    Wm = chol(Wm2);
end

% Initial guess
m0 = 0*mref*ones(mesh.nzc,1);
d0 = get1DMTfwd(m0,model,mesh,f);

% target misfit
chifact = 1;
target  = sqrt(2*chifact*nd);

%% Set up inverse problem for SPGL
% J = get1DMTJ(model,mesh,ops,f,e,1);

modelpred = model;
modelpred.m = m0;

params.Wd     = Wd;
params.Wm     = Wm; 

funFwd        = @(varargin) getDataJvSPGL(modelpred,mesh,f,varargin);
funFwdGN      = @(m)        get1DMTfwd(m,model,mesh,f);
funPenalty    = @phid;
funPenaltyH   = @(varargin) skipFirstOutput(@phid,varargin); 
funPrim       = @phim;     
funPrimH      = @(varargin) skipFirstOutput(@phim,varargin); 
funObjective  = @(m,params) objectiveTradeOff(dobs,m,@phid,@phim,funFwdGN,params);
funObjectiveH = @(varargin) skipFirstOutput(funObjective,varargin);
funDual       = @(m,Wm,params)    norm((Wm')\m,2);
project       = @WProjector;



if testit
    fprintf('Testing Forward Model \n'); 
    fun = @(m) get1DMTfwd(m,model,mesh,f);
    testGrad(fun,m0,ntest);
    fprintf('\n\n');
    
    fprintf('Testing Penalty Function \n'); 
    fun = @(r) funPenalty(r,params);
    testGrad(fun,rand(2*numel(f),1),ntest);
    fprintf('\n\n');
    
    fprintf('Testing Penalty Hessian \n'); 
    fun = @(r) funPenaltyH(r,params);
    testGrad(fun,rand(2*numel(f),1),ntest);
    fprintf('\n\n');
    
    fprintf('Testing Primary Function \n'); 
    fun = @(m) funPrim(m,Wm,params);
    testGrad(fun,rand(numel(m0),1),ntest);
    fprintf('\n\n');
    
    fprintf('Testing Primary Hessian \n');
    fun = @(m) funPrimH(m,Wm,params);
    testGrad(fun,rand(numel(m0),1),ntest);
    fprintf('\n\n');

    params.beta = 0; 
    fprintf('Testing Objective Function \n'); 
    fun = @(m) funObjective(m,params);
    testGrad(fun,rand(size(m0)),ntest);
    fprintf('\n\n');
    
    fprintf('Testing Objective Function Hessian \n'); 
    fun = @(m) funObjectiveH(m,params);
    testGrad(fun,rand(size(m0)),ntest);
    fprintf('\n\n');
    params = rmfield(params,'beta'); 
    keyboard
end

%% put options and params on structure for SPGLGeneral

% function handles
options.primal_norm = funPrim;
options.dual_norm   = funDual;
options.project     = project;
options.funPenalty  = funPenalty;
%options.nPrevVals   = 10;
% options.stepMax     = 1e1;
% weights
options.weights     = Wm; 
options.iterations  = 100; 
% options.bpTol       = 1e-4; 
%options.decTol      = 1e-3; 
options.saveOptPath = 1;
options.saveX       = 1;


%% Run Inversion Using SPGLGeneral
[mpred,r,g,info] = spgl1General(funFwd, dobs, 0, target, m0, options, params);
save('spgl1GeneralResults','mpred','r','g','info'); 

%% Gauss newton params 
% parameters for Gauss Newton with cooling
% params.beta0       = 1000; 
params.coolingrate = 4;
params.coolingfact = 1/2;


%% Run Inversion Using Gauss-Newton with Cooling on Beta
 
[mpred_gn, r_gn, info_gn] = InvGNcoolBeta(funFwdGN, dobs, m0, target,options, params);


%% Plot Results
% data
figure
semilogx(f,dtrue(1:nf),'-',f,dtrue(1:nf)+r(1:nf),'o',f,dtrue(1:nf)+r_gn(1:nf),'v','linewidth',2,'markersize',8);
hold on; 
semilogx(f,dtrue(nf+1:end),'--',f,dtrue(nf+1:end)+r(nf+1:end),'s',f,dtrue(nf+1:end)+r_gn(nf+1:end),'^','linewidth',2,'markersize',8);
xlabel('Frequency (Hz)', 'fontsize',16)
ylabel('Datum','fontsize',16)
set(gca,'fontsize',16)
legend('True, Real','Predicted SPGL, Real', 'Predicted GN, Real', 'True, Imag','Predicted SPGL, Imag','Predicted GN, Imag'); 

%% model
figure
semilogx(mesh.zc,mtrue,'-',mesh.zc,mpred+mref,'o',mesh.zc,mpred_gn+mref,'s','linewidth',2)
xlabel('Depth, z (m)', 'fontsize',16)
ylabel('m(z) + m_{ref}','fontsize',16)
set(gca,'fontsize',16)
axis([1e-1,1e6,min(mtrue)-0.5,max(mtrue)+0.5])
legend('True','Recovered SPGL','Recovered GN','location','northwest')


