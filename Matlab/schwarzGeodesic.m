% Schwarzschild geodesic without constants of motion
% Spherical standard
% Alexander Nieddu

function dx = schwarzGeodesic(x)

% Mass
M = 1;

% Initial conditions for position and velocity
% -> (t, r, theta, phi, dt, dr, dtheta, dphi)
t       = x(1);
r       = x(2);
theta   = x(3);
phi     = x(4);
dt      = x(5);
dr      = x(6);
dtheta  = x(7);
dphi    = x(8);

% Schwarzschild radius R = 2GM / c^2, G=M=c=1
R = 2;
w = 1 - R/r;
w_ = R/r^2;
v = 1/w;
v_ = -R/(R-r)^2;

% ODE derived from the geodesic equation for the schwarzschild metric 
dx(1) = dt;
dx(2) = dr;
dx(3) = dtheta;
dx(4) = dphi;
dx(5) = -1/w* w_ *dt*dr;
    
dx(6) = -(1/(2*v)*(v_)*dr^2 - r/v * dtheta^2 - (r*sin(theta)^2)/v * dphi^2 + ...
        1/(2*v)*w_*dt^2);
    
dx(7) = -(2/r*dtheta*dr - sin(theta)*cos(theta)*dphi^2);

dx(8) = -(2/r*dphi*dr + 2*cot(theta)*dphi*dtheta);

dx = dx';

end