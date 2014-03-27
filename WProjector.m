function x = twoProjector(y,Wm,tau,params)

if isempty(Wm)
    Wm = 1;
end

r = norm(Wm*y,2);

if (r <= tau)
    x = y;
    return
else
    x = y/r*tau;
end
