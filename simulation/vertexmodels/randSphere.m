function [X,Y,Z,THETA,PHI] = randSphere(N, R)

    dtheta = 0.01;
    theta = 0:dtheta:pi;
    Dist = abs(sin(theta));
    Dist = Dist./sum(Dist);
    %plot(theta, Dist)

    % inverse of cumulative distribution
    CumDist = cumsum(Dist);
    dx = 10^(-5);
    x = 0:dx:1;
    Y = interp1(CumDist, theta, x,'spline');
    %plot(x,Y)

    % dthetabin = 0.1;
    % bins = 0:dthetabin:(max(theta)+dthetabin);
    % n = histc(Y(I),bins);
    % n = n./sum(n);
    % n = n*dtheta/dthetabin;
    % 
    % bar(bins + dthetabin/2,n)
    % hold on
    % plot(theta,Dist,'-r')
    % hold off

    % draw from inverse of cumulative dist
    X = rand([N 2]);
    I = round(X(:,1)/dx) + 1;

    THETA = Y(I)'; 
    PHI = X(:,2)*2*pi;

    % convert back to cartesian
    X = R*cos(PHI).*sin(THETA);
    Y = R*sin(PHI).*sin(THETA);
    Z = R*cos(THETA);

end