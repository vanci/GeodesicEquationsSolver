function testJacobianTemplate()
    caseLocalSum();
    caseGeodesicEquationsOnSphere();
    fprintf('test JacobianTemplate\tPASS\n');
end

function caseLocalSum()
    X = cell(4,1);
    X(1:4) = {1};
    
    functions = cell(4,1);
    functions(1:3) = {@mysum};
    functions{4} = @localSum;
    
    dependency = cell(4,2);
    dependency(1:3,1) = {[1 2]; [2 3]; [3 4]};
    dependency{4,2} = 1:3;
    
    [z,JTW] = JacobianTemplate(functions, X, eye(2), dependency, []);
    assert( 0 == norm( [4;4] - z ) );
    assert( 0 == norm( [1;2;1;0] - JTW(:,1) ));
    assert( 0 == norm( [0;1;2;1] - JTW(:,2) ));
end

function y = mysum(x,Extra)
    y = sum(x);
end

function z = localSum(x,Extra)
    z = zeros(2,1);
    z = [mysum(x(1:2)); mysum(x(2:3))];
end

function caseGeodesicEquationsOnSphere()
    TOL = 1e-6;
    N = 400;
    x0 = [3*pi/2;pi/3]; xT = [pi;0];
    u0 = generateInitialValue(x0,xT,N);
 
    Extra.dim = 2; Extra.x0 = x0; Extra.xT = xT;
    tic
    [ftp, Jtp] = GeodesicEquationsOnSphereWithTemplate(u0, Extra);
    toc
    tic;
    [ffd, Jfd] = GeodesicEquationsOnSphereWithFD(u0,Extra);
    toc
    
    assert( 0 == max(ffd-ftp) )
    assert( max(max(Jfd-Jtp)) <= TOL )
end



function u0 = generateInitialValue(x0,xT,N)
    dim = size(x0,1);
    V = repmat( (xT - x0)/(N+1), N, 1 );
    steps = reshape(repmat(1:N, dim, 1), [], 1);
    u0 = repmat(x0, N, 1) + V.*steps;
end

function [f, J] = GeodesicEquationsOnSphereWithFD(u, Extra)
    M = size(u,1);
    J = zeros(M);
    F = @(x) GeodesicEquationsOnSphere(x,Extra);
    f = F(u);
    for k = 1:M
        xp = u; xp(k) = xp(k) + 1e-6;
        J(:,k) = 1e6*(F(xp) - f);
    end
end
