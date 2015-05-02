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
        C2(:,i,j) = g\C1(:,i,j);
    end
end