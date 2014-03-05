function [d, J, E, B, ops] = get1DMTfwd(mesh, model, f, sb,se,solveopt)

dz  = mesh.dz;
nzc = mesh.nzc; 
nz  = mesh.nz; 
nf  = numel(f); 

% Boundary Condidtions can be better!!!
if nargin < 6
    solveopt = 'E';
    if nargin < 5
        se = sparse(nz,1);
        if nargin < 4
            sb = sparse(nzc,1);
            sb(1) = 1;
        end
    end
end

%%
% Extract model parameters
mu  = model.mu;
sig = model.sigma;
omega = 2*pi*f;

if isfield(model,'transform')
    transform  = @(sigma) model.transform(sigma);
    itransform = @(sigma) model.itransform(sigma);
    transformderiv = @(sigma) model.transformderiv(sigma); 
else
    transform  = @(sigma) sigma;
    itransform = @(sigma) sigma;
    transformderiv = @(sigma) 1;
end


% Set up operators
G    = ddx(nzc);
Av   = ave(nzc); 
Linv = sdiag(1./dz); 
Mmu  = sdiag(dz./mu);
Msig = sdiag(Av'*(sig.*dz));
M    = sdiag(Av'*dz);

if nargout > 4
    ops.G  = G;
    ops.Av = Av;
    ops.Linv = Linv;
    ops.Mmu  = Mmu;
    ops.Msig = Msig;
    ops.M    = M;
end


% set up matrix system
switch upper(solveopt)
    case 'BOTH'
        getA = @(w) [1i*w*speye(nzc)  Linv*G; ...
             -G'*Linv'*Mmu    -Msig];
         s = @(w) [sb; M*se];
    case 'E'
         getA = @(w) (G'*Linv'*Mmu*Linv*G - 1i*w*Msig);
         s = @(w) 1i*w*M*se + G'*Linv'*Mmu*sb;
end
         

% ititialize data vector
d = zeros(nf,1);

if nargout > 1
    E = zeros(nz,nf);
    if nargout > 2
        B = zeros(nzc,nf);
    end
end

if nargout > 1
    J = zeros(nz*nf,nzc);
end

% loop over freq
for j = 1:nf
    A   = getA(omega(j));
    rhs = s(omega(j));
    u   = A\rhs;
    switch upper(solveopt)
        case 'BOTH'
            b  = u(1:nzc);
            e  = u(nzc+1:end);
        case 'E'
            e = u;
            b = (1i*omega(j)).^-1*(sb - Linv*G*e);
    end
    
    d(j)  = e(1);
    if nargout > 1
        grad = -1i*omega(j)*sdiag(e)*Av'*sdiag(dz)*sdiag(transformderiv(sig));
        Ag   = A\grad;
        J(j,:) = Ag(1,:);
        
        if nargout > 2
            E(:,j) = e;
            if nargout > 3
                B(:,j) = b;
            end
        end
    end
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