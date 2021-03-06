% function mesh = getMesh(n,sig, f)
%
% getMesh retrieves a 1D logarithmically spaced mesh with n number of 
% cells, and maximum depth = 3 skin depths, where
%
% D ~  500 / sqrt(sigma*min(f))
%
% is the skin depth
%
%
% INPUTS: 
%   n   : number of cells
%   sig : conductivity used to compute skin depth
%   f   : frequency
%
% OUTPUTS:
%   mesh : sturcture with mesh parameters including
%      dz  : node spacing
%      z   : node positions
%      zc  : cell-center positions
%      nz  : number of z-nodes
%      nzc : number of z-centers
% 
% 
% Lindsey J. Heagy
% lheagy@eos.ubc.ca 
% last modified: April 5, 2014
% ------------------------------------------------------------------------%

function mesh = getMesh(n,sig,f,zmin,zmax)

if nargin < 4
    zmin = 500/sqrt(max(sig)*max(f));
    zmin = zmin/2;
end

if nargin < 5
    zmax = 500/sqrt(min(sig)*min(f));
    zmax = 2*zmax;
end

zc = logspace(log10(zmin),log10(zmax),n+1);
dz = diff(zc);
z  = [cumsum(dz)]-dz(1);
zc = z(1:end-1) + diff(z)/2;

% put everything on a structure
mesh.dz  = dz;
mesh.z   = z;
mesh.zc  = zc;
mesh.nz  = numel(z);
mesh.nzc = numel(zc);

end