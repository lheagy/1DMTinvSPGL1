%% Test Gradient
%
% Use the Taylor series approximation to show that the addition of the
% gradient increases the approximation of a nearby point to O(2), while
% only having function information is an O(1) approximation.
%
%       bool = testGrad(fun,x);
% 
% Here, bool indicates if the test is passed, x is a random point, and fun
% is the function that has the form:
%
%       [f, g] = fun(x);
%
% Or outputs a cell with the function value in the first cell and gradient
% in the second cell.
%
% For example:
%
%       fun = @(x)({sin(x),cos(x)});
%       testGrad(fun,rand(1));
%
% 
% Rowan Cockett
% 07-Nov-2012
% University of British Columbia
% rcockett@eos.ubc.ca
function bool = testGrad(fun,x,n)
fun = @(z)(fcn(fun,z));%wrap input to ensure correct format


% Test the derivatives
[f g] = fun(x);
v = randn(size(x));

if nargin < 3
    n = 7;
end
tab = nan(1,4);
for i=1:n
    h = 10^(-i+2);
    fnew = fun(x+h*v);
    tab(i,:) = [i, norm(fnew - f), norm(fnew - f - h*g*v) nan];
    if i > 1
        tab(i,4) = log10(tab(i-1,3)./tab(i,3));
    end
    table('f(x+h*v) = f(x) + h*g(x)*v + O(3)',tab,...
        'Header',{'log(h)','O(1)?','O(2)?','Order'},...
        'PrintRow',true,'format',{'-%i','%1.1e','%1.1e','%1.4f'},'align','center')
end

num = 3;% only look at the order of last few
est_dcr = mean(tab(end-num:end,4));
tru_dcr = 1.9;%should be 2nd order
if est_dcr > tru_dcr
    fprintf('\n\nGradient Operator is working. \nMean Decrease: %4.2f\n\n',est_dcr);
else
    fprintf('\n\n%s Check the Gradient Operator. %s\n\nMean Decrease: %4.2f\n\n',repmat('*',1,10),repmat('*',1,10),est_dcr);
end

if nargout > 0
    bool = est_dcr > tru_dcr;
end

end

function [f g] = fcn(fun,x)
try
    [f, g] = fun(x);
catch e
    % maybe it is in a cell
    fg = fun(x);
    f = fg{1};
    g = fg{2};
end
end