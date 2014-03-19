function [d, ops, E, B] = get1DMTfwd(model, mesh, f, sb,se)

dz  = mesh.dz;
nzc = mesh.nzc;
nz  = mesh.nz;
nf  = numel(f);

mu0 = 4e-7*pi;

% Boundary Condidtions can be better!!!

if nargin < 5
    se = sparse(nz,1);
    if nargin < 4
        sb = sparse(nzc,1);
        sb(1) = 1;
    end
end


%%
% Extract model parameters
mu  = model.mu;
m   = model.m;
omega = 2*pi*f; 

if isfield(model,'transform')
    transform  = @(m) model.transform(m);
    itransform = @(sigma) model.itransform(sigma);
    transformderiv = @(m) model.transformderiv(m);
else
    transform  = @(m) m;
    itransform = @(sigma) sigma;
    transformderiv = @(m) 1;
end
    
ops  = getOps(model,mesh);
Linv = ops.Linv;
G    = ops.G;
Mmu  = ops.Mmu;
Msig = ops.Msig;
M    = ops.M;


s = @(w) 1i*w*M*se + G'*Linv'*Mmu*sb;


% ititialize data vector
d = zeros(nf,1);

if nargout > 1
    E = zeros(nz,nf);
    if nargout > 2
        B = zeros(nzc,nf);
    end
end

% loop over freq
for j = 1:nf
    A   = getA(omega(j),ops);
    rhs = s(omega(j));
    u   = A\rhs;
    
    e = u;
    b = (1i*omega(j)).^-1*(sb - Linv*G*e);

    d(j)  = e(1);
    
    if nargout > 2
        E(:,j) = e;
        if nargout > 3
            B(:,j) = b;
        end
    end
end

end
