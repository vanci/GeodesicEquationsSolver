function testGetChristoffelSymbols()
TOL = 1e-15;
for phi = -pi/3:pi/12:pi/3
    [C1, C2] = getChristoffelSymbols(metricSphere(phi), metricSphere_i(phi));
    if 0 ~= phi
        assert( 3 == nnz(C1) );
        assert( 3 == nnz(C2) );
        assert( sin(phi)*cos(phi) == C1(2,1,1) );
        assert( -sin(phi)*cos(phi) == C1(1,2,1) );
        assert( -sin(phi)*cos(phi) == C1(1,1,2) );
        
        assert( sin(phi)*cos(phi) == C2(2,1,1) );
        assert( norm(-tan(phi)-C2(1,2,1)) < TOL );
        assert( norm(-tan(phi)-C2(1,1,2)) < TOL );
    else
        assert( 0 == nnz(C1) );
        assert( 0 == nnz(C2) );
    end
end
fprintf('test getChristoffelSymbols\tPASS\n');
end

function g = metricSphere(phi)
    g = zeros(2,2);
    g(1,1) = cos(phi)^2;
    g(2,2) = 1;
end

function g_i = metricSphere_i(phi)
    g_i = zeros(2,2,2);
    g_i(2,1,1) = -2*sin(phi)*cos(phi);
end