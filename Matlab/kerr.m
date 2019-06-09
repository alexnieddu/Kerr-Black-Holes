% Kerr geodesic using constants of motion
% Boyer Lindquist
% Alexander Nieddu

function dx = kerr(t, x, const)
 
% Angular momentum a = J/M (0..1)
a = .7;
% Mass
M = 1;
 
% Initial parameter (r, theta, phi, t, p_r, p_theta)
r       = x(1);
theta   = x(2);
phi     = x(3);
t       = x(4);
p_r     = x(5);
p_theta = x(6);

% Constants of motion
mu = const(1);
E = const(2);
L = const(3);
k = const(4);
 
% Abbreviations
Sigma = r^2 + a^2 * cos(theta)^2;
Delta = r^2 - 2 * M * r + a^2;

% ODE from Fuerst & Wu
dx(1) = (p_r * Delta) / Sigma;
dx(2) = p_theta / Sigma;
dx(3) = (2*a*r*E + (Sigma - 2*r)*L*csc(theta)^2) / (Sigma * Delta);
dx(4) = E + (2*r*(r^2 + a^2)*E - 2*a*r*L) / (Sigma * Delta);
dx(5) = 1/(Sigma * Delta) * (((r^2 + a^2)*mu - k)*(r-1) + r*Delta*mu + ...
        2*r*(r^2 + a^2)*E^2 - 2*a*E*L) - (2*p_r^2*(r-1)) / (Sigma);
dx(6) = (sin(theta)*cos(theta)) / (Sigma) * ((L^2)/(sin(theta)^4) - a^2*(E^2 + mu));
 
dx = dx';
 
end