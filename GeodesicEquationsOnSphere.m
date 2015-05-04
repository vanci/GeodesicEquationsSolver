function z = GeodesicEquationsOnSphere(x, Extra)
    x0 = Extra.x0; xT = Extra.xT;
    M = size(x,1); dim = size(x0,1); N = M/dim;
    z = zeros(M,1);
    pre = 1:dim;  curr = pre + dim;  pos = curr + dim;
    
    X = [x0;x;xT];
    for k = 1:N
        vl = X(curr) - X(pre);
        vm = (X(pos) - X(pre))/2;
        vr = X(pos) - X(curr);
        z(pre) = vr - vl + getCrossTerm(X(curr), vm);
        pre = curr; curr = pos; pos = pos + dim;
    end
end

function a = getCrossTerm(x, v)
    dim = size(x,1);
    a = zeros(dim,1);
    [~, C2] = getChristoffelSymbols(metricSphere(x), metricSphere_i(x));
    
    for k = 1:dim
        a(k) = v'*C2(:,:,k)'*v;
    end
end

function [C1, C2] = getChristoffelSymbols(g, g_i)
    dim = size(g,1);

    % compute Christoffel symbol of the first kind
    C1 = zeros(dim,dim,dim);
    Idx_gamma = zeros(dim^3,1);
    Idx_g1 = zeros(dim^3,1);
    Idx_g2 = zeros(dim^3,1);
    Idx_g3 = zeros(dim^3,1);
    for k = 1:dim
        for i = 1:dim
            for j = 1:dim
                idx = k+dim*i+dim^2*j;
                Idx_gamma(idx) = idx;
                Idx_g1(idx) = idx;
                Idx_g2(idx) = k+dim*j+dim^2*i;
                Idx_g3(idx) = i+dim*j+dim^2*k;
                C1(k,i,j) =  (g_i(i,k,j) + g_i(j,i,k) - g_i(k,i,j)) / 2;
            end
        end
    end

    % compute Christoffel symbol of the second kind
    C2 = zeros(dim,dim,dim);
    for i = 1:dim
        for j = 1:dim
            C2(i,j,:) = g\C1(:,i,j);
        end
    end
end

function g = metricSphere(u)
    g = zeros(2,2);
    g(1,1) = cos(u(2))^2;
    g(2,2) = 1;
end

function g_i = metricSphere_i(u)
    g_i = zeros(2,2,2);
    g_i(2,1,1) = -2*sin(u(2))*cos(u(2));
end