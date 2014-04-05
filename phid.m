% function  [f,g] = phid(r,params)
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

function  [f,g] = phid(r,params)

Wd = params.Wd;
wr = Wd*r;

w2r = Wd'*wr;

f = norm(wr,2);
g = w2r/f;

