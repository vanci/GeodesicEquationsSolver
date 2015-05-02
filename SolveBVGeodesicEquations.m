function [X, V] = SolveBVGeodesicEquations( metric, metric_i, x0, xT, N)
    options = optimoptions('fsolve','TolFun',1e-8);

    %% Finite Difference
    u0 = generateInitialValue(x0,xT,N);
    F = @(u)finiteDifferenceEquations(metric, metric_i, x0, u, xT);
    u = fsolve(F,u0,options);
    X = reshape(u,size(x0,1),[]);
    V = 0;
end

function u0 = generateInitialValue(x0,xT,N)
    dim = size(x0,1);
    V = repmat( (xT - x0)/(N+1), N, 1 );
    steps = reshape(repmat(1:N, dim, 1), [], 1);
    u0 = repmat(x0, N, 1) + V.*steps;
end

function z = finiteDifferenceEquations(metric, metric_i, x0, x, xT)
    M = size(x,1); dim = size(x0,1); N = M/dim;
    z = zeros(M,1);
    pre = 1:dim;  curr = pre + dim;  pos = curr + dim;
    
    X = [x0;x;xT];
    for k = 1:N
        vl = X(curr) - X(pre);
        vm = (X(pos) - X(pre))/2;
        vr = X(pos) - X(curr);
        z(pre) = vr - vl + getCrossTerm(metric, metric_i, X(curr), vm);
        pre = curr; curr = pos; pos = pos + dim;
    end
end

function a = getCrossTerm(metric, metric_i, x, v)
    dim = size(x,1);
    a = zeros(dim,1);
    [~, C2] = getChristoffelSymbols(metric(x), metric_i(x));
    
    for k = 1:dim
        a(k) = v'*reshape(C2(k,:,:),dim,dim)'*v;
    end
end
