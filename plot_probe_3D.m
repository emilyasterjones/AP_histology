function plot_probe_3D(tv, probes)

cmap = colormap('lines');

% Plot probe trajectories
figure('Name','Probe trajectories');
axes_atlas = axes;
plotBrainGrid([],axes_atlas);
set(axes_atlas,'ZDir','reverse');
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])
h = rotate3d(gca);
h.Enable = 'on';

% plot probe trajectories
for p = 1:length(probes)
    line(probes{p}.trajectory_coords([1,end],1),...
        probes{p}.trajectory_coords([1,end],3),...
        probes{p}.trajectory_coords([1,end],2), ...
        'color', cmap(p,:), 'linewidth',2)
end

end