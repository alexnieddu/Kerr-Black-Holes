% Tranformation to cartesian coordinates
% Alexander Nieddu

function x_ = cartesian(x, a)

l = length(x(:,1));
x_ = ones(l,3);

r = x(:,1);
theta = x(:,2);
phi = x(:,3);

x_(:,1) = (sqrt(r.^2 + a^2)).*sin(theta).*cos(phi);
x_(:,2) = (sqrt(r.^2 + a^2)).*sin(theta).*sin(phi);
x_(:,3) = r.*cos(theta);

%x_ = x_';

end