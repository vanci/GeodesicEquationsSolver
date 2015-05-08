function [X, V] = SolveBVGeodesicEquationsOnSphere( x0, xT, N)
    TolFun = 1e-8;

    %% Finite Difference
    u0 = generateInitialValue(x0,xT,N);
    Extra.x0 = x0; Extra.xT = xT; Extra.N = N;
    
    options = optimoptions('fsolve','TolFun',TolFun,'Jacobian', 'off');
    F = @(x) GeodesicEquationsOnSphere(x,Extra);
    x = fsolve(F,u0,options);
    
    X = reshape(x,size(x0,1),[]);
    V = 0;
end

function u0 = generateInitialValue(x0,xT,N)
    dim = size(x0,1);
    V = repmat( (xT - x0)/(N+1), N, 1 );
    steps = reshape(repmat(1:N, dim, 1), [], 1);
    u0 = repmat(x0, N, 1) + V.*steps;
end