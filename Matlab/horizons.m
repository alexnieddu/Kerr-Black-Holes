% Plotting Horizons
% Plots inner, outer horizon and ergosphere in section
% Alexander Nieddu

a=0:0.01:1;
F(length(a)) = struct('cdata',[],'colormap',[]);
for j = 1:length(a)+1

    a = j/(102);
    M = 1;
    theta = 0:0.1:2*pi;

    % Outer ergosphere
    r = M + sqrt(M.^2 - a.^2*cos(theta).^2);
    polarplot(theta, r)
    hold on;

    % Inner ergosphere
    r = M - sqrt(M.^2 - a.^2*cos(theta).^2);
    polarplot(theta, r)
    hold on;

    % Outer horizon
    r = M + sqrt(M.^2 - a.^2*ones(length(theta)));
    polarplot(theta, r)
    hold on;

    % Inner horizon
    r = M - sqrt(M.^2 - a.^2*ones(length(theta)));
    polarplot(theta, r)
    hold on;

    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaAxisUnits = 'radians';
    ax.ThetaZeroLocation = 'top';
    ax.NextPlot = 'replaceChildren';
    F(j) = getframe;
end

v = VideoWriter('out.avi');
v.FrameRate = 60;
open(v);
writeVideo(v, F);