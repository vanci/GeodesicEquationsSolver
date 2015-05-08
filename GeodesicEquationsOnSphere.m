function z = GeodesicEquationsOnSphere(x, Extra)
    x0 = Extra.x0; xT = Extra.xT;
    M = size(x,1); dim = Extra.dim; N = M/dim;
    z = zeros(M,1);
    pre = 1:dim;  curr = pre + dim;  pos = curr + dim;
    
    X = [x0;x;xT];
    for n = 1:N
        u = X(curr);
        vl = X(curr) - X(pre);
        vm = (X(pos) - X(pre))/2;
        vr = X(pos) - X(curr);
        
        %% compute nonzeros of metric tensor
        gInv11 = cos(u(2))^-2;
        nz_idx = 1+2^2;
        g_i_nz = -2*sin(u(2))*cos(u(2));

        %% compute acceleration
        r = acceleration(vm,dim,nz_idx,g_i_nz);
        
        z(pre) = vr - vl + [gInv11*r(1); r(2)];
        pre = curr; curr = pos; pos = pos + dim;
    end
end

function r = acceleration(vm,dim,nz_idx,g_i_nz)
    r = 0*vm;
    for i = 1:dim
        for j = 1:dim
            C1Column = 0*vm;
            for k = 1:dim
                if nz_idx == (k+dim*(i-1)+dim^2*(j-1))
                    C1Column(k) = g_i_nz/2;
                elseif nz_idx == (k+dim*(j-1)+dim^2*(i-1))
                    C1Column(k) = g_i_nz/2;
                elseif nz_idx == (i+dim*(j-1)+dim^2*(k-1))
                    C1Column(k) = -g_i_nz/2;
                end
            end
            r = r + C1Column*vm(i)*vm(j);
        end
    end
end