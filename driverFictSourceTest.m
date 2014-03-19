% fictSourceTest
% 
% clear all
% close all
%clc

f = 1 ; %logspace(-2,5,10);
N = 11;

sigma = @(x) x.^2+1;
mu    = @(x) cos(x) + 2;

b = @(x) x.*(x-1).*mu(x); % note b(0) = b(1) = 0;
e = @(x) cos(2*pi*x);
% The system
% 1i*w*b + e’ = s1
% (1/mu * b)’ - sigma*e = s2
% fictitious source
s1 = @(x) 1i*2*pi*f*b(x) - 2*pi*sin(2*pi*x);
s2 = @(x) 2*x-1 - sigma(x).*e(x);

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
    h = rand(n,1); L = sum(h);
    h = h/L;
    mesh = getMesh(n,0,h);
    z = mesh.z;
    zc = mesh.zc; 
    model.sigma = sigma(zc);
    model.mu    = mu(zc);
    
    sb = s1(zc);
    se = s2(z);
    
    [d, ops, E, B] = get1DMTfwd(model,mesh,f,sb,se,'E');
    
    errB2(i) = norm(ops.Linv\(B - b(zc)));
    errBi(i) = norm(B - b(zc),'inf');
    errE2(i) = norm(ops.Linv\(1/2*(E(1:n)+E(2:n+1)) - e(zc)));
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