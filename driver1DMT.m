% 1DFDEM Inversion

clear all
close all
clc

addpath('../spgl1');

%% Mesh Parameters
nc   = 5000; % number of core cells
np   = 20; % number of padding cells
dz   = 1; 
mu0  = 4e-7*pi;
mesh = getMesh(nc,np,dz);
%% Model

sig = [1e-2 1e-4 1e-1 1e-2];
m   = log(sig); 
d   = [0  100 300 800]; 

sref = 1e-4;

model = getLayerModel(mesh,d,m);
model.transform = @(m) exp(m);
model.itransform = @(sigma) log(sigma);
model.transformderiv = @(m) sdiag(exp(m));
model.mref = model.itransform(sref);

%% Analytical Solutions
k1 = @(sig,f) sqrt(-1i*mu0*sig*2*pi*f);
Bana = @(sig,f,z) exp(1i*k1(sig,f)*-z);
Eana = @(sig,f,z) -2*pi*f/k1(sig,f)*exp(1i*k1(sig,f)*-z);

%% forward model
se = sparse(mesh.nz,1);
sb = sparse(mesh.nzc,1);
sb(1) = 1;
f = logspace(0,5,30);
[dobs,ops,e] = get1DMTfwd(model,mesh,f,sb,se);

%% Set up weights for objective functions

% Data weight
perc = 0.03;
flr  = 0.01*mean(abs(dobs(:)));
Wd   = perc*abs(dobs(:)) + flr;
Wd   = sdiag(1./Wd);

%% Set up inverse problem for SPGL
% J = get1DMTJ(model,mesh,ops,f,e);

funFwd = @(varargin) getDataJvSPGL(model,mesh,f,varargin);

%[mpred,r,g,info] = spgl1General(funFwd, dobs, ~, nf, m0, options, params)