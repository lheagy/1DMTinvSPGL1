% 1DFDEM Inversion

clear all
close all
clc

addpath('../spgl1');

%% Mesh Parameters
nc = 100; % number of core cells
np = 50; % number of padding cells
dz = 1; 

mesh = getMesh(nc,np,dz);

%% Model

sig = [1e-2 1e-2 1e-2];
d   = [0    100   200]; 

model = getLayerModel(mesh,d,sig);

%% forward model

f = logspace(-2,2,10);
[dobs,J] = get1DMTfwd(mesh,model,f);
dana = anaMT1Dsolu([d' sig'], f);