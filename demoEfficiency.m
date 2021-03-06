function demoEfficiency()
%% this script compare efficiency of different ways of computing Jacobian
    caseNewtonMethod(2000);
    %for N = 1000:500:5000
    %    caseSphereJacobian(N);
    %end
end

function caseNewtonMethod(N)
    x0 = [pi/2;pi/3]; xT = [3*pi/2;-pi/3];
    dim = length(x0);

    Extra.x0 = x0; Extra.xT = xT; Extra.N = N; Extra.dim = dim;
    Extra.JPIs = [];
    u0 = generateInitialValue(x0,xT,N);
    
    tic
    F = @(x) GeodesicEquationsOnSphere(x,Extra);
    [u1, r1] = fsolve(F,u0);
    tokeiFSolve = toc
    
    r2 = r1;
    u2 = u0;
    tic
    while norm(r2) >= norm(r1)
        [r2, J, Extra.JPI] = GeodesicEquationsOnSphereWithSparseAD(u2, Extra);
        u2 = u2 - J\r2;
    end
    tokeiSADSolve = toc
end

%% caseSphereJacobian compute Jacobian matrix of geodesic equations on sphere
%  N is the number of nodes
function caseSphereJacobian(N)
    x0 = [pi/2;pi/3]; xT = [3*pi/2;-pi/3];
    dim = length(x0);
    
    maxIter = 3;
    diff = zeros(maxIter,5);
    res = zeros(maxIter,1);
    tokeiNT = zeros(1, maxIter);
    tokei = zeros(6,maxIter);
    Extra.x0 = x0; Extra.xT = xT; Extra.N = N; Extra.dim = dim;
    Extra.JPIs = [];
    u = generateInitialValue(x0,xT,N);
    
    for k = 1:maxIter
        
        % finite difference
        profile -memory on
        [ffd, Jfd] = GeodesicEquationsOnSphereWithFD(u,Extra);
        profreport
        profsave(profile('info'),sprintf('SphereWithFD\\N_%d_iter_%d',N,k))
        tic
        [ffd, Jfd] = GeodesicEquationsOnSphereWithFD(u,Extra);
        tokei(1,k) = toc;

        % forward mode
        profile -memory on
        [ffm, Jfm] = GeodesicEquationsOnSphereWithForwardModeAD(u,Extra);
        profreport
        profsave(profile('info'),sprintf('SphereWithFM\\N_%d_iter_%d',N,k))
        tic
        [ffm, Jfm] = GeodesicEquationsOnSphereWithForwardModeAD(u,Extra);
        tokei(2,k) = toc;
        
        % reverse mode
%         profile -memory on
%         [frm, Jrm] = GeodesicEquationsOnSphereWithReverseModeAD(u,Extra);
%         profreport
%         profsave(profile('info'),sprintf('SphereWithRM\\N_%d_iter_%d',N,k))
%         tic
%         [frm, Jrm] = GeodesicEquationsOnSphereWithReverseModeAD(u,Extra);
%         tokei(3,k) = toc;
        
        % bi-coloring sparse mode
        profile -memory on
        [fsp, Jsp, Extra.JPI] = GeodesicEquationsOnSphereWithSparseAD(u, Extra);
        profreport
        profsave(profile('info'),sprintf('SphereWithSAD\\N_%d_iter_%d',N,k))
        tic
        [fsp, Jsp, Extra.JPI] = GeodesicEquationsOnSphereWithSparseAD(u, Extra);
        tokei(4,k) = toc;
        
        % template mode
        profile -memory on
        [ftp, Jtp] = GeodesicEquationsOnSphereWithTemplate(u, Extra);
        profreport
        profsave(profile('info'),sprintf('SphereWithTP\\N_%d_iter_%d',N,k))
        tic
        [ftp, Jtp] = GeodesicEquationsOnSphereWithTemplate(u, Extra);
        tokei(5,k) = toc;
        
        % sparse template mode
        profile -memory on
        [fstp, Jstp] = GeodesicEquationsOnSphereWithSparseTemplate(u, Extra);
        profreport
        profsave(profile('info'),sprintf('SphereWithSTP\\N_%d_iter_%d',N,k))
        tic
        [fstp, Jstp] = GeodesicEquationsOnSphereWithSparseTemplate(u, Extra);
        tokei(6,k) = toc;
        
        
        diff(k,1) = max(max(Jfd-Jfm));
%        diff(k,2) = max(max(Jfd-Jrm));
        diff(k,3) = max(max(Jfd-Jsp));
        diff(k,4) = max(max(Jfd-Jtp));
        diff(k,5) = max(max(Jfd-Jstp));
        res(k) = norm(ffd);
        assert( 0 == norm(ffd-ffm) )
%        assert( 0 == norm(ffd-frm) )
        assert( 0 == norm(ffd-fsp) )
        assert( 0 == norm(ffd-ftp) )
        assert( 0 == norm(ffd-fstp) )
        % Newton step
        tic;
        u = u - Jfd\ffd;
        tokeiNT(k) = toc;
    end
    
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
    
    filename = sprintf('tokei_N_%d.mat',N);
    save( filename, 'tokei', 'tokeiNT' );
    tokei
    tokeiNT
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