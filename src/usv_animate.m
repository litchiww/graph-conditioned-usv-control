function usv_animate(R, fig_dir, opts)

if nargin < 3, opts = struct; end
if ~isfield(opts,'fps'),        opts.fps = 20; end
if ~isfield(opts,'duration_s'), opts.duration_s = []; end
if ~isfield(opts,'gif'),        opts.gif = true; end
if ~isfield(opts,'mp4'),        opts.mp4 = true; end

cfg = usv_config(); layout = cfg.scenario_layout();
USV_COLORS = [0.122 0.467 0.706; 0.839 0.153 0.157; 0.173 0.627 0.173];

a0 = R.C0_baseline.arr;
a2 = R.C2_full.arr;
T  = min(a0.t(end), a2.t(end));
if isempty(opts.duration_s), opts.duration_s = T; end
n_frames = round(opts.fps * opts.duration_s);
ts = linspace(0, T, n_frames);

fig = figure('Visible','off','Position',[100 100 1500 720],'Color','w');
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

axB = nexttile(tl); local_setup_axes(axB, layout, 'C0 baseline');
axF = nexttile(tl); local_setup_axes(axF, layout, 'C2 full method (ours)');

hL_B = gobjects(3,1); hL_F = gobjects(3,1);
hT_B = gobjects(3,1); hT_F = gobjects(3,1);
for j = 1:3
    hL_B(j) = animatedline(axB,'Color',USV_COLORS(j,:),'LineWidth',1.8);
    hL_F(j) = animatedline(axF,'Color',USV_COLORS(j,:),'LineWidth',1.8);
    hT_B(j) = patch(axB, NaN, NaN, USV_COLORS(j,:),'EdgeColor','k','LineWidth',0.8);
    hT_F(j) = patch(axF, NaN, NaN, USV_COLORS(j,:),'EdgeColor','k','LineWidth',0.8);
end

hText = annotation(fig,'textbox',[0.45 0.95 0.1 0.04], 'String','t = 0.0 s', ...
    'FontSize',13,'FontWeight','bold','HorizontalAlignment','center', ...
    'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5]);

if opts.mp4
    mp4_path = fullfile(fig_dir, 'animation_baseline_vs_ours.mp4');
    vw = VideoWriter(mp4_path, 'MPEG-4');
    vw.FrameRate = opts.fps; vw.Quality = 92;
    open(vw);
end
if opts.gif
    gif_path = fullfile(fig_dir, 'animation_baseline_vs_ours.gif');
end

for fi = 1:n_frames
    tnow = ts(fi);
    iB = find(a0.t <= tnow, 1, 'last');
    iF = find(a2.t <= tnow, 1, 'last');
    for j = 1:3

        [Xb, Yb, psib] = local_phys(a0.t(1:iB), a0.x1(1:iB,j), a0.x2(1:iB,j), j, layout);
        clearpoints(hL_B(j)); addpoints(hL_B(j), Xb, Yb);
        [xs, ys] = local_tri(Xb(end), Yb(end), psib(end), 70);
        set(hT_B(j), 'XData', xs, 'YData', ys);

        [Xf, Yf, psif] = local_phys(a2.t(1:iF), a2.x1(1:iF,j), a2.x2(1:iF,j), j, layout);
        clearpoints(hL_F(j)); addpoints(hL_F(j), Xf, Yf);
        [xs, ys] = local_tri(Xf(end), Yf(end), psif(end), 70);
        set(hT_F(j), 'XData', xs, 'YData', ys);
    end
    set(hText, 'String', sprintf('t = %.1f s', tnow));
    drawnow limitrate;
    frame = getframe(fig);
    if opts.mp4, writeVideo(vw, frame); end
    if opts.gif
        [imind, cm] = rgb2ind(frame2im(frame), 256);
        if fi == 1
            imwrite(imind, cm, gif_path, 'gif','LoopCount', Inf,'DelayTime', 1/opts.fps);
        else
            imwrite(imind, cm, gif_path, 'gif','WriteMode','append','DelayTime', 1/opts.fps);
        end
    end
    if mod(fi, max(1,round(n_frames/10))) == 0
        fprintf('  frame %d / %d\n', fi, n_frames);
    end
end
if opts.mp4, close(vw); fprintf('  wrote %s\n', mp4_path); end
if opts.gif, fprintf('  wrote %s\n', gif_path); end
close(fig);
end

function local_setup_axes(ax, layout, label)
hold(ax,'on'); ax.Color = [0.92 0.95 0.98];
for r = 1:size(layout.turbines,1)
    xc = layout.turbines(r,1); yc = layout.turbines(r,2);
    rectangle(ax,'Position',[xc-layout.safety_radius, yc-layout.safety_radius, ...
                             2*layout.safety_radius, 2*layout.safety_radius], ...
              'Curvature',[1 1],'EdgeColor',[0.5 0.5 0.5],'LineStyle',':','LineWidth',0.7);
    plot(ax, xc, yc, 'o', 'Color',[0.27 0.27 0.27],'MarkerFaceColor',[0.27 0.27 0.27],'MarkerSize',6);
    plot(ax, [xc xc],[yc-25 yc+25],'Color',[0.27 0.27 0.27],'LineWidth',1.2);
    plot(ax, [xc-35 xc+35],[yc yc], 'Color',[0.4 0.4 0.4],'LineWidth',1.0);
end
for yc = layout.lane_centers_y
    plot(ax,[0 4200],[yc yc],'--','Color',[0.6 0.6 0.6],'LineWidth',0.6);
end
xlabel(ax,'East (m)'); ylabel(ax,'North (m)');
title(ax, label, 'FontSize',12);
xlim(ax, layout.x_corridor); ylim(ax, layout.y_corridor);
axis(ax,'equal'); grid(ax,'on'); ax.GridLineStyle=':'; ax.GridAlpha=0.3;
end

function [X,Y,psi] = local_phys(t,x1,x2,row,layout)
Y = layout.lane_centers_y(row) + x1*layout.length_scale;
X = layout.path_progress_rate * t;
psi = atan2(x2, 1.0);
end

function [xs, ys] = local_tri(x, y, psi, sz)
ang0 = psi + [pi/2, pi/2 + 2*pi/3, pi/2 + 4*pi/3];
xs = x + sz*cos(ang0); ys = y + sz*sin(ang0);
end
