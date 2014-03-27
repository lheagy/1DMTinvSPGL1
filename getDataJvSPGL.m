% function [out] = getDataJvSPGL(model,mesh,f,varargin)
%
% Wrapper for get1DMTfwd to be consistent with SPGLGeneral

function [out] = getDataJvSPGL(model,mesh,f,varargin)

m = varargin{1}{1};
v  = varargin{1}{2};

if ~isempty(v)
    [~,J] = get1DMTfwd(m,model,mesh,f);
    out = J'*v;
else
    [dobs] = get1DMTfwd(m,model,mesh,f);
    out = dobs;
end
