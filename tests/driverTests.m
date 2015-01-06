% fictSourceTest
% 
clear all
close all
clc

f = 100; %logspace(-2,5,10);
N = 7;

mu0 = 4e-7*pi;
w = 2*pi*f;

model.transform = @(m) exp(m);
model.itransform = @(m) log(m);

r = mu0*w; 

% analytical functions
% fields
u1a = @(x) (x-1).^3;
u2a = @(x) cos(pi*x);
u3a = @(x) (x-1).^2;
u4a = @(x) sin(pi*x);

m   = @(x) ones(size(x)); %sin(2*pi*x);

s1a = @(x) 2*(x-1).^2;
s2a = @(x) -(1+pi)*sin(pi*x);
s3a = @(x) 2*(x-1) + r*exp(m(x)).*cos(pi*x);
s4a = @(x) pi*cos(pi*x) + r*exp(m(x)).*(x-1).^3;

e = @(x) u1a(x) + 1i*u2a(x);
b = @(x) u3a(x) + 1i*u4a(x); 
%%
fprintf('log2(n)    errB2   oB2     errBi   oBi     errE2   oE2     errEi   oEi     \n');   

errB2 = zeros(N-2,1);
errBi = zeros(N-2,1);
errE2 = zeros(N-2,1);
errEi = zeros(N-2,1);

% model.transform  = @(sigma) log(sigma);
% model.itransform = @(sigma) exp(sigma);
% model.transformderiv = @(sigma) diag(sparse(exp(sigma)));

for i = 1:N-2
    n = 2^(3+i);
    h = ones(n,1); L = sum(h);
    h = h/L;
    mesh = getMesh(n,0,h);
    z = mesh.z;
    zc = mesh.zc;
    
    mod = m(zc);
    model.mu = mu0; 
    nz = mesh.nz;
    
    s1 = s1a(z);
    s2 = s2a(z);
    s3 = s3a(z);
    s4 = s4a(z);
    
    q = zeros(nz*4,1);
    q(1:4:end) = s1;
    q(2:4:end) = s2;
    q(3:4:end) = s3;
    q(4:4:end) = s4;
    
    d = get1DMTfwd(mod, model, mesh, f, [], q);
    E = d(1:4:end)+1i*d(2:4:end);
    B = d(3:4:end)+1i*d(4:4:end);
    
    errB2(i) = norm(sdiag(mesh.dz)*(1/2*(B(1:n)+B(2:n+1)) - b(zc)));
    errBi(i) = norm(B - b(z),'inf');
    errE2(i) = norm(sdiag(mesh.dz)*(1/2*(E(1:n)+E(2:n+1)) - e(zc)));
    errEi(i) = norm(E - e(z),'inf');
    
    if i < 2
        oB2 = NaN;
        oBi = NaN;
        oE2 = NaN;
        oEi = NaN;
    else
        oB2 = log2(errB2(i-1)/errB2(i));
        oBi = log2(errBi(i-1)/errBi(i));
        oE2 = log2(errE2(i-1)/errE2(i));
        oEi = log2(errEi(i-1)/errEi(i));
    end
    
    fprintf('%5i   %1.2e  %2.2f  %1.2e  %2.2f  %1.2e  %2.2f  %1.2e  %2.2f  \n', log2(n), errB2(i), oB2, errBi(i), oBi, errE2(i), oE2, errEi(i), oEi);
end