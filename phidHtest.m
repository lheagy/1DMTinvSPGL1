% Wrapper so I can test Hessian
function [g,H] = phidHtest(r,params)
    [~,g,H] = phid(r,params);
end