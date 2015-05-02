function testJacobianTemplate()
    caseLocalSum();
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
    fprintf('test JacobianTemplate\tPASS\n');
end

function y = mysum(x,Extra)
    y = sum(x);
end

function z = localSum(x,Extra)
    z = zeros(2,1);
    z = [mysum(x(1:2)); mysum(x(2:3))];
end

function c2 = getChristoffelSymbols(x, Extra)
end

function g = metricSphere(u,Extra)
    g = zeros(2,2);
    g(1,1) = cos(u(2))^2;
    g(2,2) = 1;
end

function g_i = metricSphere_i(u,Extra)
    g_i = zeros(2,2,2);
    g_i(2,1,1) = -2*sin(u(2))*cos(u(2));
end

