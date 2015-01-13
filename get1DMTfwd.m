% function [dobs,J] = get1DMTfwd(m,model,mesh,f,v,q)
%
% get1DMTfwd solves the 1D magnetotelluric problem
%
% Follows: 
%   E. Haber, U. Ascher, and D. Oldenburg, On optimization techniques for
%   solving nonlinear inverse problems, Inverse Problems 16 (2000)
%   1263-1280
%
% INPUTS:
%   deltam : m - mref
%   model  : structure containing fields 
%               mref : reference model (optional, default = 0) 
%               transform : model transform that maps the model to
%                               electrical conductivity
%               itransform : inverse transform - maps from conductivity to
%                               model space
%               transformderiv : derivative of the model transform that
%                               maps from model to electrical conductivity
% 
%   mesh   : structure containing fields
%               dz  : cell spacing
%               z   : nodal mesh
%               zc  : cell centered mesh
%               nz  : number of nodes
%               nzc : number of cell centers
%
%   f      : vector of frequencies 
%   v      : (optional) vector to multiply (J*v) 
%   q      : (optional) source term, default = [1,0,0,...,0]'
%
% OUTPUTS:
%   dobs   : observed data, [er;ei]
%   J      : sensitivity matrix 
%
% Lindsey J. Heagy
% lheagy@eos.ubc.ca
%


function [dobs,J] = get1DMTfwd(deltam,model,mesh,f,v,q)

if isfield(model,'mref')
    mref = model.mref;
else
    mref = 0;
end

if ~isfield(model,'transform')
    model.transform = @(m) m     ; 
    model.transformderiv = @(m) 1; 
end


mod = deltam+mref;

A = getApaper(mod,mesh,f);

nz = mesh.nz;
nf = numel(f);

if nargin < 5
    q = sparse(4*nz*nf,1);
    q([1:4*nz:end]) = 1;
end

u = A\q;

er = u(1:4:end);
ei = u(2:4:end);
erobs = er(1:nz:end);
eiobs = ei(1:nz:end);
dobs = [erobs; eiobs];

if nargout > 1
    if nargin > 4
        J = getJ(mod,mesh,f,u,A,model,v);
    else
        J = getJ(mod,mesh,f,u,A,model);
    end
    Jr = J(1:4:end,:); Jr = Jr(1:nz:end,:);
    Ji = J(2:4:end,:); Ji = Ji(1:nz:end,:);
    J = [Jr;Ji];
end


end

function A = getApaper(m,mesh,f)

nf = numel(f);

mu0 = 4e-7*pi;

w = 2*pi*f;
r = mu0*w;

zc = mesh.zc;
nzc = mesh.nzc;
nz  = mesh.nz;

B0 = zeros(2,4);
B0(1,3) = 1;
B0(2,4) = 1;

B1= zeros(2,4);
B1(1,1) = 1;
B1(2,2) = 1;

A = [];

for j = 1:nf
    Rblks = B0;
    Sblks = [];
    for i = 1:nzc
        %ind = 4*(i-1)+1:4*i;
        R     = getR(zc(i),r(j),m(i));
        Rblks = blkdiag(Rblks,R);
        S     = getS(zc(i),r(j),m(i));
        Sblks = blkdiag(Sblks,S);
    end
    Sblks = blkdiag(Sblks,B1);
    Sblks = [sparse(2,4*nz); Sblks];
    Rblks = [Rblks; sparse(2,4*nz)];
    
    Af = Sblks + Rblks;
    A = blkdiag(A,Af);
end
    
end

function G = getG(m,mesh,f,u,model)

nf = numel(f);

mu0 = 4e-7*pi;

w = 2*pi*f;
r = mu0*w;

dz = mesh.dz; 
nzc = mesh.nzc;
nz  = mesh.nz;

u1 = u(1:4:end);
u2 = u(2:4:end);

G = [];
for j = 1:nf
    Gk = sparse(2,nzc);
    u1f = u1((j-1)*nz+1:j*nz);
    u2f = u2((j-1)*nz+1:j*nz);
    for i = 1:nzc
        gki = getGki(dz(i),r(j),m(i),u1f(i:i+1),u2f(i:i+1),model);
        Gki = zeros(4,nzc);
        Gki(:,i) = gki; 
        Gk  = [Gk;Gki]; 
    end
    Gk = [Gk;sparse(2,nzc)];
    G  = [G;Gk];
end

end


function J = getJ(m,mesh,f,u,A,model,v) %deltam,mesh,f,u,A,model
G = getG(m,mesh,f,u,model);
if nargin>6
    G = G*v;
end
J = -A\G;

end

function T = getT(r,m)
rem = r*exp(m);
T = zeros(4,4);
T(1,3) = 1; T(2,4) = 1;
T(3,2) = -rem;
T(4,1) = rem;
end

function S = getS(h,r,m)
I = speye(4);
T = getT(r,m);
S = -h.^-1*I-T;
end

function R = getR(h,r,m)
I = speye(4);
T = getT(r,m);
R = h.^-1*I-T;
end

function Gki = getGki(~,r,m,u1,u2,model)
transformderiv = model.transformderiv;
Gki = sparse(4,1);
Gki(3) = r*transformderiv(m)*sum(u2)/2;
Gki(4) = -r*transformderiv(m)*sum(u1)/2;
end
