% function V = sdiag(v)
%
% Lindsey J. Heagy - credit: Eldad Haber
%
% March 26, 2014
function V = sdiag(v)
V = diag(sparse(v));
end