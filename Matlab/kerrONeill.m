% Kerr geodesic using constants of motion from O Neill
% Boyer Lindquist
% Alexander Nieddu

function dx = kerrONeill(t, x, const)
 
% Angular momentum a = J/M (0..1)
a = .84;
% Mass
M = 1;
 
% Initial parameter (r, theta, phi, t)
r       = x(1);
theta   = x(2);
phi     = x(3);
t       = x(4);

% Constants of motion
mu = const(1);
E = const(2);
L = const(3);
Q = const(4);
 
% Abbreviations
Sigma = r^2 + a^2 * cos(theta)^2;
Delta = r^2 - 2 * M * r + a^2;
K = Q + (L - a*E)^2;
X = a^2*(E^2 + mu) - L^2 - Q;
R = (E^2 + mu)*r^4 - 2*M*mu*r^3 + X*r^2 + 2*M*K*r -a^2*Q;
O = Q + cos(theta)^2*(a^2*(E^2 + mu) - L^2/sin(theta)^2);

% ODE from O Neill Chapter 4.2.2 Geodesics
dx(1) = -0*sqrt(R / Sigma^4);
dx(2) = 0*sqrt(O / Sigma^4);
dx(3) = ((L - a*sin(theta)^2*E) / sin(theta)^2 + a*((r^2 + a^2)*E - L*a)/(Delta) ) / Sigma^2;
dx(4) = (a*(L - a*E*sin(theta)^2) + (r^2 + a^2)*((r^2 + a^2)*E - L*a)/(Delta)) / Sigma^2;
 
dx = dx';
 
end