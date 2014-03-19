function model = getLayerModel(mesh,d,sig,mu)

nd  = numel(d);
zc  = mesh.zc;
nzc = mesh.nzc;

if nargin < 4;
    mu = 4e-7*pi*ones(nd,1);
end

if abs(d(1)) > 1e-6
    warning('getLayerModel:TopLayerHeight','Top Layer not at z = 0'); 
end


msigma = sig(end)*ones(nzc,1);
mmu    = mu(end)*ones(nzc,1);

for i = nd-1:-1:1
    ic = zc <= d(i+1);
    msigma(ic) = sig(i);
    mmu(ic)    = mu(i);
end

model.m = msigma;
model.mu    = mmu;
