function [d] = get1DMTfwd(h, f, model, display)

if nargin < 3
    display = false;
end

if display; fprintf('1DMT FORWARD MODEL \n'); end
%%
% mesh parameters
n = numel(h);
nf = numel(f); 

% Extract model parameters
mu  = model.mu;
sig = model.sigma;

omega = 2*pi*f;

% compute skin depths 
skdpth = sqrt(500./([min(sig) max(sig)].*[min(f) max(f)]));
if display; fprintf(' min skdepth = %f, max skdepth = %f \n', min(skdpth), max(skdpth));end

% Set up operators
G    = ddx(n);
Av   = ave(n); 
Linv = sdiag(1./h);
Mmu  = sdiag(h./mu);
Msig = sdiag(Av'*(sig.*h));
M    = sdiag(Av'*h);

% set up matrix system
getA = @(w) [1i*w*speye(n)  Linv*G; ...
             -G'*Linv'*Mmu  -Msig ];

% boundary conditions - this can be better!
bb    = [zeros(n,1); zeros(n+1,1)];
bb(1) = 1;

% ititialize data vector
d = zeros(nf,1);

% loop over freq
for j = 1:nf
    A = getA(omega(j));
    be = A\bb;
    b  = be(1:n);
    e  = be(n+1:end);
    d(j)  = e(1);
end

end

function G = ddx(n)
G = spdiags(ones(n+1,1)*[-1 1], [0,1],n,n+1);
end

function V = sdiag(v)
V = diag(sparse(v));
end

function Av = ave(n)
Av = spdiags(ones(n+1,1)*[1/2 1/2],[0 1],n,n+1);
end