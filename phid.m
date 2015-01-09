% function  [f,g,H] = phid(r,params)
% 
% data misfit function and derivative
%
%   f = || Wd r ||_2
%
%   g = (Wd' Wd r) / f
%
% INPUTS:
%   r : residual
%   params : structure with field Wd containing data weighting matrix
%
% OUTPUTS:
%   f : data misfit
%   g : gradient of data misfit
%
% Lindsey J. Heagy
% last modified: March 26, 2014

function  [f,g,H] = phid(r,params)

Wd = params.Wd;

if isfield(params,'eps')
    eps = params.eps;
else
    eps = 1e-10;
end


Wdr = Wd*r; 
f   = norm(Wdr,2);

if nargout > 1
    if f < eps
        g = zeros(numel(r),1);
    else
        Wd2r = Wd'*Wdr;
        g = Wd2r/f;
    end
    
    if nargout > 2 
        if f < eps
            H = Wd'*Wd; % this is cheating!
        else
            WdtWd = Wd'*Wd;
            rrt   = r*r';
            H     = -(WdtWd*rrt*WdtWd)/f^3 + WdtWd/f ;
        end
    end
end


