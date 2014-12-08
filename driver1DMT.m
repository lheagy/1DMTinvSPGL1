% 1DFDEM Inversion
% 
% Driver file for the 1D magnetotelluric forward and inverse problems.
%
%
% Lindsey J. Heagy
% lheagy@eos.ubc.ca
% last modified: April 5, 2014
%
% ------------------------------------------------------------------------%

%% House-keeping 
clear all
close all
clc

% add path to SPGL1general
addpath(genpath('../spgl1'));

%% Parameters

% mesh
nc   = 90;      % number of padding cells

% frequency 
frange = [0 2]; % log10 of fmin, fmax
nf   = 10;   % number of frequencies 
nd   = nf*2; % number of data

% reference model
mref = -2;      % log of reference sigma
mu0  = 4e-7*pi; % magnetic permeability

%% Mesh
mesh = getMesh(nc,exp(mref),10.^frange);

%% Model
eta = 4;
t = log(mesh.zc);

mb = -2;
mtrue = mb+2*exp(-(t-0.025).^2/2)-1*exp(-(t-3.5).^2/3);
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
noise = percn*randn(nd,1).*abs(dtrue);
dobs = dtrue+noise;

% Data weight
perc = 0.03;
flr  = 0.003;
flrr  = flr*mean(abs(dobs(1:nf)));
flri  = flr*mean(abs(dobs(nf+1:end)));
Wd   = perc*abs(dobs) + [flrr*ones(nf,1);flri*ones(nf,1)];
Wd   = sdiag(1./Wd);

% Model Weights
alphas = 1;
alphaz = 2;
alpha2z = 8; 
Wms = speye(mesh.nzc);

t   = log(mesh.z);
dt  = log(mesh.dz);
dts = (dt(1:end-1) + dt(2:end))/2;
Wmd = sdiag([1./dts])*spdiags([-1*ones(nzc,1), ones(nzc,1)],-1:0,nzc,nzc);
Wmd(1,1) = 0;

% sdiag([1./dts]).^2*
Wmd2 = spdiags([-1*ones(nzc,1), 2*ones(nzc,1), -1*ones(nzc,1)],[-1:1],nzc,nzc);
Wmd2(1,[1:2]) = 0; 
Wmd2(end,[end-1:end]) = 0; 

Wm  = [alphas*Wms; alphaz*Wmd; alpha2z*Wmd2];


[n,m] = size(Wm);
if n~=m
    Wm2 = Wm'*Wm;
    Wm = chol(Wm2);
end

% Initial guess
m0 = mref*ones(mesh.nzc,1);
d0 = get1DMTfwd(m0,model,mesh,f);

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

% function handles
options.primal_norm = funPrim;
options.dual_norm   = funDual;
options.project     = project;
options.funPenalty  = funPenalty;
%options.nPrevVals   = 10;
% options.stepMax     = 1e1;
% weights
options.weights     = Wm; 
options.iterations  = 50; 
% options.optTol      = target*0.01;
options.decTol      = 1; 
params.Wd           = Wd;



%% Run Inversion Using SPGLGeneral
[mpred,r,g,info] = spgl1General(funFwd, dobs, 0, target, m0, options, params);


%% Plot Results
% data
figure
semilogx(f,dtrue(1:nf),f,dtrue(1:nf)+r(1:nf),'o','linewidth',2,'markersize',8); hold on; semilogx(f,dtrue(nf+1:end),'--',f,dtrue(nf+1:end)+r(nf+1:end),'s','linewidth',2,'markersize',8);
xlabel('Frequency (Hz)', 'fontsize',16)
ylabel('Datum','fontsize',16)
set(gca,'fontsize',16)
legend('True Data, Real Component','Predicted Data, Real Component','True Data, Complex Component','Predicted Data, Complex Component')


% model
figure
semilogx(mesh.zc,mtrue,mesh.zc,mpred+mref,'o','linewidth',2)
xlabel('Depth, z (m)', 'fontsize',16)
ylabel('m(z) + m_{ref}','fontsize',16)
set(gca,'fontsize',16)
axis([mesh.z(1) mesh.z(end),-3,0])
legend('True Model','Recovered Model','location','northwest')
