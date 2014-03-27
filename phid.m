% function  [f,g] = phid(r,params)
% data misfit function
% Lindsey J. Heagy
% last modified: March 26, 2014

function  [f,g] = phid(r,params)

Wd = params.Wd;
wr = Wd*r;

w2r = Wd'*wr;

f = norm(wr,2);
g = w2r/f;

