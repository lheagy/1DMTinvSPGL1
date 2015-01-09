function [f,g,H] = phim(m,Wm,params)

if isfield(params,'eps')
    eps = params.eps;
else
    eps = 1e-10;
end

Wmm = Wm*m;
f   = norm(Wmm,2);

if nargout > 1

    if f < eps
        g = zeros(numel(m),1);
    else
        g = Wm'*(Wm*m)/f;
    end
    if nargout > 2

        if f < eps
            H = Wm'*Wm; % This is cheating!
            %sparse(numel(m));
        else
            WmtWm = Wm'*Wm;
            mmt   = m*m';
            H     = -(WmtWm*mmt*WmtWm)/f^3 + WmtWm/f ;
        end
    end
end

end
