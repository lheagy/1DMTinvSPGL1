function mesh = getMesh(nc,np,h)
if nargin < 3
    h = 1;
end

dz  = h.*[ones(nc,1); 1.3.^[1:np]'];
z  = [0; cumsum(dz)];
zc = z(1:end-1) + diff(z)/2;

mesh.dz  = dz;
mesh.z   = z; 
mesh.zc  = zc;
mesh.nz  = numel(z);
mesh.nzc = numel(zc);

end