% Wrapper so I can test Hessian
function [g,H] = skipFirstOutput(fun,varargin)

if nargin < 2
    [~,g,H] = fun(inputs{1});
else
    inputs = varargin{:};
    n = numel(inputs);
    if n == 2
        [~,g,H] = fun(inputs{1},inputs{2});
    elseif n == 3
        [~,g,H] = fun(inputs{1},inputs{2},inputs{3});
    else
        [~,g,H] = fun(varargin);
    end
end