function testSolveBVGeodesicEquationsOnSphere()
    caseSphere();
end

function caseSphere()
    TOL = 1e-4;
    h = 2*pi/100;
    N = ceil(2*pi/h);

    xb0 = [3*pi/2;pi/3]; xbT = [pi;0];
    [Xb, Vb] = SolveBVGeodesicEquationsOnSphere(xb0, xbT, N);
    
    xc0 = [pi;0]; xcT = [pi/2;-pi/3];
    [Xc, Vc] = SolveBVGeodesicEquationsOnSphere(xc0, xcT, N);
    
    xa0 = [pi;0]; xaT = [0;0];
    [Xa, Va] = SolveBVGeodesicEquationsOnSphere(xa0, xaT, N);
    
    close all
    sphere; hold on;
    verifySphere(Xa, 'r', TOL);
    verifySphere(Xb, 'g', TOL);
    verifySphere(Xc, 'b', TOL); hold off;
    
    diff = zeros(16,1);
    Extra.x0 = xb0; Extra.xT = xbT; Extra.N = N;
    u = generateInitialValue(xb0,xbT,N);
    for k = 1:16
        [ffd, Jfd] = GeodesicEquationsOnSphereWithFD(u,Extra);
        [fad, Jad] = GeodesicEquationsOnSphereWithAD(u,Extra);
        diff(k) = max(max(Jfd-Jad));
        res(k) = norm(fad);
        assert( 0 == norm(ffd-fad) )
        u = u - Jfd\ffd;
    end
    figure; subplot(2,1,1)
    plot(diff)
    subplot(2,1,2)
    plot(res);
end

function verifySphere(X, color, TOL)
    x = cos(X(1,:)).*cos(X(2,:));
    y = sin(X(1,:)).*cos(X(2,:));
    z = sin(X(2,:));
    plot3(x,y,z,color,'LineWidth',2);
    Y = [x;y;z];
    N = size(Y,2);
    % any three points on the great circle should lie in the same plane
    % with the origin
    U = Y(:,[1,ceil(N/5),ceil(3*N/4)]);
    Z = cross(U(:,1),U(:,2));
    p = dot(Z,U(:,3));
    assert( abs(p) < TOL, num2str(abs(p)));
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

function [f, J] = GeodesicEquationsOnSphereWithAD(u, Extra)
    M = size(u,1);
    AD_fun = ADfun(@GeodesicEquationsOnSphere, M);
    options = setopt('revprod', eye(M));
    [f, J] = feval(AD_fun, u, Extra, options);
end

function u0 = generateInitialValue(x0,xT,N)
    dim = size(x0,1);
    V = repmat( (xT - x0)/(N+1), N, 1 );
    steps = reshape(repmat(1:N, dim, 1), [], 1);
    u0 = repmat(x0, N, 1) + V.*steps;
end