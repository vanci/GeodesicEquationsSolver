function [f, J, JPIs] = GeodesicEquationsOnSphereWithSparseTemplate(u, Extra)
    dim = Extra.dim; N = Extra.N;
    
    functions = cell(N+1,1);
    functions{1} = @geodesicEquationLeftBoundary;
    functions(2:(end-2)) = {@geodesicEquation};
    functions{end-1} = @geodesicEquationRightBoundary;
    functions{end} = @combineResult;
    
    dependency = cell(N+1,2);
    dependency{1,1} = [1 2];
    d = repmat( [1;2;3], N-2, 1) + reshape( repmat(0:(N-3), 3, 1), [], 1);
    dependency(2:(end-2),1) = num2cell(reshape(d,3,[])',2);
    dependency{end-1,1} = [N-1 N];
    dependency{end,2} = 1:N;
    
    U = num2cell(reshape(u,dim,[]), 1)';
    
    [f,JTW, JPIs] = JacobianTemplateWithColoring(functions, U, speye(dim*N), Extra.JPIs, dependency, Extra);
    J = JTW';
end

function z = geodesicEquationLeftBoundary(x, Extra)
    dim = Extra.dim;
    z = zeros(dim,1);
    curr = 1:dim; pos = curr + dim;
    vl = x(curr) - Extra.x0;
    vm = (x(pos) - Extra.x0)/2;
    vr = x(pos) - x(curr);
    
    %% compute nonzeros of metric tensor
    gInv11 = cos(x(curr(2)))^-2;
    nz_idx = 1+2^2;
    g_i_nz = -2*sin(x(curr(2)))*cos(x(curr(2)));
    
    %% compute acceleration
    r = acceleration(vm,dim,nz_idx,g_i_nz);
    
    z = vr - vl + [gInv11*r(1); r(2)];
end

function z = geodesicEquation(x, Extra)
    dim = Extra.dim;
    z = zeros(dim,1);
    pre = 1:dim; curr = pre + dim; pos = curr + dim;
    vl = x(curr) - x(pre);
    vm = (x(pos) - x(pre))/2;
    vr = x(pos) - x(curr);
    
    %% compute nonzeros of metric tensor
    gInv11 = cos(x(curr(2)))^-2;
    nz_idx = 1+2^2;
    g_i_nz = -2*sin(x(curr(2)))*cos(x(curr(2)));
    
    %% compute acceleration
    r = acceleration(vm,dim,nz_idx,g_i_nz);
    
    z = vr - vl + [gInv11*r(1); r(2)];
end

function z = geodesicEquationRightBoundary(x, Extra)
    dim = Extra.dim;
    z = zeros(Extra.dim,1);
    pre = 1:dim; curr = pre + dim;
    vl = x(curr) - x(pre);
    vm = (Extra.xT - x(pre))/2;
    vr = Extra.xT - x(curr);
    
    %% compute nonzeros of metric tensor
    gInv11 = cos(x(curr(2)))^-2;
    nz_idx = 1+2^2;
    g_i_nz = -2*sin(x(curr(2)))*cos(x(curr(2)));
    
    %% compute acceleration
    r = acceleration(vm,dim,nz_idx,g_i_nz);
    
    z = vr - vl + [gInv11*r(1); r(2)];
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

function z = combineResult(x, Extra)
    z = x;
end