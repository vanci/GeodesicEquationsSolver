function testSolveBVGeodesicEquationsOnSphere()
    %caseSphereSolution();
    caseSphereJacobian();
end

function caseSphereJacobian()
    N = 100;
    
    xb0 = [pi/2;pi/3]; xbT = [3*pi/2;-pi/3];
    dim = length(xb0);
    
    maxIter = 3;
    diff = zeros(maxIter,5);
    res = zeros(maxIter,1);
    tokeiNT = zeros(1, maxIter);
    tokei = zeros(6,maxIter);
    Extra.x0 = xb0; Extra.xT = xbT; Extra.N = N; Extra.dim = dim;
    Extra.JPIs = [];
    u = generateInitialValue(xb0,xbT,N);
    
    for k = 1:maxIter
        
        % finite difference
        tic
        [ffd, Jfd] = GeodesicEquationsOnSphereWithFD(u,Extra);
        tokei(1,k) = toc;
        
        % forward mode
        tic
        [ffm, Jfm] = GeodesicEquationsOnSphereWithForwardModeAD(u,Extra);
        tokei(2,k) = toc;
        
        % reverse mode
        tic
        [frm, Jrm] = GeodesicEquationsOnSphereWithReverseModeAD(u,Extra);
        tokei(3,k) = toc;
        
        % bi-coloring sparse mode
        tic
        [fsp, Jsp, Extra.JPI] = GeodesicEquationsOnSphereWithSparseAD(u, Extra);
        tokei(4,k) = toc;
        
        % template mode
        tic
        [ftp, Jtp] = GeodesicEquationsOnSphereWithTemplate(u, Extra);
        tokei(5,k) = toc;
        
        % sparse template mode
        tic
        [fstp, Jstp] = GeodesicEquationsOnSphereWithSparseTemplate(u, Extra);
        tokei(6,k) = toc;
        
        
        diff(k,1) = max(max(Jfd-Jfm));
        diff(k,2) = max(max(Jfd-Jrm));
        diff(k,3) = max(max(Jfd-Jsp));
        diff(k,4) = max(max(Jfd-Jtp));
        diff(k,5) = max(max(Jfd-Jstp));
        res(k) = norm(ffd);
        assert( 0 == norm(ffd-ffm) )
        assert( 0 == norm(ffd-frm) )
        assert( 0 == norm(ffd-fsp) )
        assert( 0 == norm(ffd-ftp) )
        assert( 0 == norm(ffd-fstp) )
        % Newton step
        tic;
        u = u - Jfd\ffd;
        tokeiNT(k) = toc;
    end
    sphere; hold on;
    verifySphere(reshape(u,2,[]), 'r', Inf); hold off
    
    figure; subplot(6,1,1)
    plot(diff(:,1));  title('difference between FD and Forward Mode')
    subplot(6,1,2)
    plot(diff(:,2));  title('difference between FD and Reverse Mode')
    subplot(6,1,3);
    plot(diff(:,3));  title('difference between FD and Bi-Coloring')
    subplot(6,1,4);
    plot(diff(:,4));  title('difference between FD and Template Mode')
    subplot(6,1,5);
    plot(diff(:,5));  title('difference between FD and Sparse Template Mode')
    subplot(6,1,6);
    plot(res);        title('residual v.s. num. of newton steps');
    
    tokei
    tokeiNT
end

function caseSphereSolution()
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

function [f, J] = GeodesicEquationsOnSphereWithReverseModeAD(u, Extra)
    M = size(u,1);
    AD_fun = ADfun(@GeodesicEquationsOnSphere, M);
    options = setopt('revprod', eye(M));
    [f, J] = feval(AD_fun, u, Extra, options);
    J = J';
end

function [f, J] = GeodesicEquationsOnSphereWithForwardModeAD(u, Extra)
    M = size(u,1);
    AD_fun = ADfun(@GeodesicEquationsOnSphere, M);
    options = setopt('forwprod', eye(M));
    [f, J] = feval(AD_fun, u, Extra, options);
end

function [f, J, JPI] = GeodesicEquationsOnSphereWithSparseAD(u, Extra)
    M = size(u,1);
    if( isempty(Extra) || ~isfield(Extra,'JPI') )
        [JPI, ~] = getjpi(@GeodesicEquationsOnSphere, M, M, Extra, 'd');
    else
        JPI = Extra.JPI;
    end
    [f, J] = evalj(@GeodesicEquationsOnSphere, u, Extra, M, JPI);
end

function u0 = generateInitialValue(x0,xT,N)
    dim = size(x0,1);
    V = repmat( (xT - x0)/(N+1), N, 1 );
    steps = reshape(repmat(1:N, dim, 1), [], 1);
    u0 = repmat(x0, N, 1) + V.*steps;
end