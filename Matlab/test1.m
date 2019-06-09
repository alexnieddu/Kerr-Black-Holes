% Working simulation

a = .95; M = 1;

r0 = 6.5;
theta0 = pi/2;
phi0 =  0;

r_in = M - sqrt(M^2 - a^2);
r_out = M + sqrt(M^2 - a^2);

Delta = r0^2 - 2 * M * r0 + a^2;

% FALLING OBSERVER
% mu = -1;
% E = sqrt(2);
% L_range = 2*M*E*r_out/a - 2*M*E*r_in/a;
% L = 2*M*E*r_in/a + L_range*.4;
% L = a*E;
% Q = 0;

% PHYS.COMPUTER
mu = -1;
E = .7;
L = (2*M*a*E + sqrt(r0*Delta*((E^2 - 1)*r0 + 2*M))) / (2*M - r0);
Q = 0;

% STEIN
% mu = 0;
% E = 1;
% L = -(r0^3 - 3*M*r0^2 + a^2*r0 + a^2*M)/(a*(r0 - M));
% Q = -r0^3*(r0^3 - 6*M*r0^2 + 9*M^2*r0 + 4*a^2*M)/(a^2*(r0 - M)^2);

mu = -1;
E = 0.956545;
L = -0.830327;
Q = 13.4126;

const = [mu, E, L, Q];

x0 = [ r0 theta0 phi0 0];

% Integration
[t, res] = ode45(@(t, x) kerrONeill(t, x, const), [0 15000], x0);
%[t, res] = runge4(@(t, x) kerrONeill(t, x, const), [0 10], x0, 1e-6);
cart = cartesian(res(:,1:3),a);


% Plot geodesic
ax = 6;
plot3(cart(:,1),cart(:,2),cart(:,3));
title(["Kerr geodesic"]);
dim = [.1 .7 .3 .2];
str = {"Constants of motion:","E=" num2str(E) ", L= " num2str(L)};
annotation('textbox',dim,'String',str,'FitBoxToText','on')
xlabel("x");
ylabel("y");
zlabel("z");
hold on;
% Plot testparticle
plot3(cart(end,1),cart(end,2),cart(end,3), '-r.', 'MarkerSize', 10)
hold on;
% Plot black hole as dot with sphere
plot3(0,0,0, '-k.', 'MarkerSize',20)
[x y z] = sphere;
h = surfl(r_in*x, r_in*y, r_in*z);
shading interp
set(h,'FaceColor',[1 0 0], 'FaceAlpha', 0.85)
h = surfl(r_out*x, r_out*y, r_out*z); 
set(h,'FaceColor',[0 1 0], 'FaceAlpha', 0.1)
shading interp
axis equal;
% axis([-ax ax -ax ax -ax ax]);