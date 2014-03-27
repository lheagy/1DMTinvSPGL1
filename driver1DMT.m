% 1DFDEM Inversion
%
% Lindsey J. Heagy
% last modified: March 26, 2014
%
%% 

clear all
close all
clc

addpath(genpath('../spgl1')); % Depends on SPGL1general

%% Mesh Parameters
nc   = 0; % number of core cells
np   = 90; % number of padding cells
dz   = 1e-1;
frange = [0 2]; % log10 of fmin, fmax

mu0  = 4e-7*pi;
mesh = getMesh(nc,np,dz);
nf   = 30;

nd   = nf*2; % number of data
%% Model

mref = -1;
eta = 4;
t = log(mesh.zc);

mtrue = -1+2*exp(-(t-2).^2/2)+-0*exp(-(t-4).^2/2);
model.mu = mu0*ones(mesh.nzc,1);

model.transform = @(m) exp(m);
model.itransform = @(sigma) log(sigma);
model.transformderiv = @(m) sdiag(exp(m));
model.mref = mref;

m = mtrue-mref;
model.m = m;

%% forward model

f = logspace(frange(1),frange(2),nf);
dtrue = get1DMTfwd(m,model,mesh,f);

%% Set up weights for objective functions

nzc = mesh.nzc; 

% add noise to data
percn = 0.03;
noise = percn*randn(nd,1).*dtrue;
dobs = dtrue+0*noise;

% Data weight
perc = 0.03;
flr  = 0.01*min(abs(dobs));
Wd   = perc*abs(dobs(:)) + flr;
Wd   = sdiag(1./Wd);

% Model Weights
alphas = 1e-4;
alphaz = 1;
Wms = alphas*speye(mesh.nzc);

t   = log(mesh.z);
dt  = diff(t);
dts = (dt(1:end-1) + dt(2:end))/2;
Wmd = sdiag([0; 1./dts])*alphaz*spdiags([-1*ones(nzc,1), ones(nzc,1)],-1:0,nzc,nzc);
Wm = [Wms]; %;Wmd];

[n,m] = size(Wm);
if n~=m
    Wm2 = Wm'*Wm;
    Wm = chol(Wm2);
end

% Initial guess
m0 = mref*ones(mesh.nzc,1);
d0 = get1DMTpaper(m0,model,mesh,f);

% target misfit
chifact = 1;
target  = chifact*nd;

%% Set up inverse problem for SPGL
% J = get1DMTJ(model,mesh,ops,f,e,1);

modelpred = model;
modelpred.m = m0;

funFwd     = @(varargin) getDataJvSPGL(modelpred,mesh,f,varargin);
funPenalty = @phid;
funPrim    = @(m,Wm,params)    norm(Wm*m,2);
funDual    = @(m,Wm,params)    norm((Wm')\m,2);
project    = @WProjector;

%% put options and params on structure

options.primal_norm = funPrim;
options.dual_norm   = funDual;
options.project     = project;
options.funPenalty  = funPenalty;
options.weights     = Wm;  
params.Wd           = Wd;

%% Feed to SPGLGeneral
[mpred,r,g,info] = spgl1General(funFwd, dobs, 0, target, m0, options, params);