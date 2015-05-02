function testSolveBVGeodesicEquations()
    caseSphere();
    % caseDeSitter();
end

function caseDeSitter()
    TOL = Inf; %1e-14;
    h = 2*pi/30;
    N = ceil(2*pi/h);
    
    deSitterSurface(-pi*0.6,pi*0.6,-pi,pi,pi/20); hold on;
    x0 = [-pi/2;0]; xT = [pi/2;0];
    [Xa, Va] = SolveBVGeodesicEquations( @metricSphere, @metricSphere_i, x0, xT, N);
    verifyDeSitter(Xa,'g',TOL);
    
    tau = pi/2; c = 0;
    x0 = [-tau; 2*atan(exp(-tau))+c]; xT = [tau; 2*atan(exp(tau))+c];
    [Xb, Vb] = SolveBVGeodesicEquations( @metricSphere, @metricSphere_i, x0, xT, N);
    verifyDeSitter(Xb,'g',TOL);
    
    hold off   
end

function caseSphere()
    TOL = 1e-6;
    h = 2*pi/100;
    N = ceil(2*pi/h);

    x0 = [pi;0]; xT = [0;0];
    [Xa, Va] = SolveBVGeodesicEquations( @metricSphere, @metricSphere_i, x0, xT, N);
    
    x0 = [pi;0]; xT = [3*pi/2;pi/3];
    [Xb, Vb] = SolveBVGeodesicEquations( @metricSphere, @metricSphere_i, x0, xT, N);
    
    x0 = [pi;0]; xT = [pi/2;-pi/3];
    [Xc, Vc] = SolveBVGeodesicEquations( @metricSphere, @metricSphere_i, x0, xT, N);
    
    
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

function verifyDeSitter(X, color, TOL)
    x = sinh(X(1,:));
    y = cosh(X(1,:)).*cos(X(2,:));
    z = cosh(X(1,:)).*sin(X(2,:));
    plot3(x,y,z,color,'LineWidth',2);
    N = size(x,2); a = ceil(N/5); b = ceil(N/2); c = ceil(2*N/3);
    d1 = [x(c)-x(b); y(c)-y(b); z(c)-z(b)];
    d2 = [x(b)-x(a); y(b)-y(a); z(b)-z(a)];
    p = dot(d1,d2)/norm(d1)/norm(d2);
    assert( abs(1-abs(p)) < TOL );
end

function deSitterSurface(lb1,ub1,lb2,ub2,h)
    [U1, U2] = meshgrid(lb1:h:ub1,lb2:h:ub2);
    X = sinh(U1);
    Y = cosh(U1).*cos(U2);
    Z = cosh(U1).*sin(U2);
    mesh(X,Y,Z);
    colormap(parula(5))
    alpha(.6)
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

function g = metricDeSitter(u)
    g = zeros(2,2);
    g(1,1) = 1;
    g(2,2) = -cosh(u(1))^2;
end

function g_i = metricDeSitter_i(u)
    g_i = zeros(2,2,2);
    g_i(1,2,2) = -2*sinh(u(1))*cosh(u(1));
end