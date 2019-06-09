% Runge Kutta Method (4. Order)
% Alexander Nieddu
function [t, y ] = runge4 ( f, tRange, x0, h )
t(1) = tRange(1);
numSteps = ( tRange(2) - tRange(1) ) / h;
y(:,1) = x0;
for k = 1 : numSteps
    k1 = f(t, y(:,k) );
    k2 = f(t, y(:,k) + k1 * h/2 );
    k3 = f(t, y(:,k) + k2 * h/2 );
    k4 = f(t, y(:,k) + k3 * h );
    t(1,k+1) = t(1,k) + h;
    y(:,k+1) = y(:,k) + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4);
end
t = t';
y = y';