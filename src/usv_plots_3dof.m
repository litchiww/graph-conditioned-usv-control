function usv_plots_3dof(hydro, fig_dir)

S = usv_style(); S.apply();
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

arr = hydro.arr;
out = hydro.out;
hp = hydro.hp;
tools = usv_metrics();

fig = S.figure('double',0.72);
tl = tiledlayout(fig,3,6,'Padding','compact','TileSpacing','compact');

ax = nexttile(tl,1,[2 2]); hold(ax,'on');
for j = 1:hp.N
    plot(ax, arr.Xref(:,j), arr.Yref(:,j), ':', 'Color', S.USV(j,:), ...
        'LineWidth', 0.9, 'HandleVisibility','off');
    plot(ax, arr.X(:,j), arr.Y(:,j), S.USVLineStyle{j}, ...
        'Color', S.USV(j,:), ...
        'LineWidth', 1.25, 'DisplayName', sprintf('USV %d',j));
end
xlabel(ax,'East position (m)'); ylabel(ax,'North position (m)');
title(ax,'(a) 3-DOF trajectory tracking','FontWeight','normal');
legend(ax,'show','Location','best','FontSize',6.6);
local_grid(ax);

ax = nexttile(tl,3,[2 2]); hold(ax,'on');
err = abs(arr.lat_err);
for j = 1:hp.N
    plot(ax, arr.t, tools.moving_rms(err(:,j), arr.t, 2.0), S.USVLineStyle{j}, ...
        'Color', S.USV(j,:), 'LineWidth', 1.1, 'HandleVisibility','off');
end
href = yline(ax, out.rms_lateral_m, 'k--');
href.HandleVisibility = 'off';
text(ax,0.97,0.90,sprintf('second-half RMS %.2f m',out.rms_lateral_m), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',7.0,'BackgroundColor','w','Margin',1.5);
xlabel(ax,'time (s)'); ylabel(ax,'2-s moving RMS error (m)');
title(ax,'(b) Physical tracking error','FontWeight','normal');
xlim(ax,[0 arr.t(end)]);
local_grid(ax);

ax = nexttile(tl,5,[2 2]); hold(ax,'on');
heading_err_deg = rad2deg(atan2(sin(arr.psiref-arr.psi), ...
    cos(arr.psiref-arr.psi)));
for j = 1:hp.N
    plot(ax, arr.t, tools.moving_rms(heading_err_deg(:,j),arr.t,2.0), S.USVLineStyle{j}, ...
        'Color', S.USV(j,:), ...
        'LineWidth', 1.1, 'HandleVisibility','off');
end
href = yline(ax,out.rms_heading_deg,'k--');
href.HandleVisibility = 'off';
text(ax,0.97,0.90,sprintf('second-half RMS %.2f deg',out.rms_heading_deg), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',7.0,'BackgroundColor','w','Margin',1.5);
xlabel(ax,'time (s)'); ylabel(ax,'2-s moving RMS (deg)');
title(ax,'(c) Heading-reference error','FontWeight','normal');
xlim(ax,[0 arr.t(end)]);
local_grid(ax);

ax = nexttile(tl,13,[1 2]); hold(ax,'on');
for j = 1:hp.N
    plot(ax,arr.t,arr.delta_deg(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:), ...
        'LineWidth',1.0,'HandleVisibility','off');
end
yline(ax,hp.delta_max_deg,'k--','HandleVisibility','off');
yline(ax,-hp.delta_max_deg,'k--','HandleVisibility','off');
ylabel(ax,'rudder equivalent $\delta$ (deg)');
xlabel(ax,'time (s)');
title(ax,'(d) Limited steering command','FontWeight','normal');
xlim(ax,[0 arr.t(end)]);
ylim(ax,1.08*hp.delta_max_deg*[-1 1]);
local_grid(ax);

ax = nexttile(tl,15,[1 2]); hold(ax,'on');
power_kw = sum(arr.power_W,2)/1000;
nwin = max(3,round(5/median(diff(arr.t))));
plot(ax,arr.t,movmean(power_kw,nwin),'Color',S.safee,'LineWidth',1.1);
yline(ax,out.mean_total_power_W/1000,'k--','HandleVisibility','off');
text(ax,0.97,0.88,sprintf('mean %.2f kW',out.mean_total_power_W/1000), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',6.8,'BackgroundColor','w','Margin',1.2);
xlabel(ax,'time (s)'); ylabel(ax,'5-s mean fleet power (kW)');
title(ax,'(e) Force/moment power estimate','FontWeight','normal');
xlim(ax,[0 arr.t(end)]); local_grid(ax);

ax = nexttile(tl,17,[1 2]); hold(ax,'on');
sep_margin = arr.min_sep - hp.safety_sep;
plot(ax,arr.t,sep_margin,'Color',[0.20 0.20 0.20],'LineWidth',1.05);
yline(ax,0,'r--','HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'separation margin (m)');
title(ax,'(f) Fleet clearance','FontWeight','normal');
xlim(ax,[0 arr.t(end)]); local_grid(ax);

exportgraphics(fig, fullfile(fig_dir,'fig_hydro_validation.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig, fullfile(fig_dir,'fig_hydro_validation.png'), 'Resolution',300);
close(fig);
fprintf('  fig_hydro_validation\n');
end

function local_grid(ax)
grid(ax,'on'); ax.GridLineStyle=':'; ax.GridAlpha=0.3;
end
