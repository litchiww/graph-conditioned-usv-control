function usv_plots_3dof_closedloop(hydro,fig_dir)

S=usv_style(); S.apply();
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
arr=hydro.arr;
out=hydro.out;
hp=hydro.hp;
tools=usv_metrics();

fig=S.figure('double',0.60);
tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

ax=nexttile(tl); hold(ax,'on');
for j=1:hp.N
    plot(ax,arr.t,tools.moving_rms(arr.lat_err(:,j),arr.t,5.0), S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.1, ...
        'DisplayName',sprintf('USV %d',j));
end
yline(ax,out.rms_lateral_m,'k--','HandleVisibility','off');
text(ax,0.62,0.88,sprintf('second-half RMS %.2f m',out.rms_lateral_m), ...
    'Units','normalized','HorizontalAlignment','center', ...
    'VerticalAlignment','top','FontSize',6.8, ...
    'BackgroundColor','w','Margin',1.5);
xlabel(ax,'physical time (s)');
ylabel(ax,'5-s moving RMS (m)');
title(ax,'(a) Closed-loop lane regulation','FontWeight','normal');
legend(ax,'show','Location','northeast','NumColumns',1,'FontSize',6.4);
xlim(ax,[0 arr.t(end)]); local_grid(ax);

ax=nexttile(tl); hold(ax,'on');
for j=1:hp.N
    plot(ax,arr.t,abs(arr.ratio(:,j)),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:), ...
        'LineWidth',1.05,'HandleVisibility','off');
end
yline(ax,1,'k--','HandleVisibility','off');
xlabel(ax,'physical time (s)');
ylabel(ax,'$|x_2^{(j)}|/F_2^{(j)}$');
title(ax,'(b) Operational-envelope ratio','FontWeight','normal');
ylim(ax,[0 1.05]); xlim(ax,[0 arr.t(end)]); local_grid(ax);
text(ax,0.97,0.93,sprintf('maximum %.3f',out.max_envelope_ratio), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',7.0,'BackgroundColor','w');

ax=nexttile(tl); hold(ax,'on');
for j=1:hp.N
    plot(ax,arr.t,arr.tauN(:,j)./hp.Nmax(j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:), ...
        'LineWidth',1.0,'HandleVisibility','off');
end
yline(ax,1,'k--','HandleVisibility','off');
yline(ax,-1,'k--','HandleVisibility','off');
xlabel(ax,'physical time (s)');
ylabel(ax,'$\tau_N^{(j)}/\tau_{N,\max}^{(j)}$');
title(ax,'(c) Direct yaw-moment command','FontWeight','normal');
lim=max(0.55,1.08*out.max_yaw_moment_ratio);
ylim(ax,[-lim lim]); xlim(ax,[0 arr.t(end)]); local_grid(ax);
text(ax,0.97,0.93,sprintf('saturation %d/%d', ...
    out.yaw_saturation_count,out.control_channel_samples), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',7.0,'BackgroundColor','w');

ax=nexttile(tl); hold(ax,'on');
show=arr.t<=2.0;
for j=1:hp.N
    plot(ax,arr.t(show),arr.eta2(show,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:), ...
        'LineWidth',1.1,'HandleVisibility','off');
end
yline(ax,hp.direction_confidence2_min,'k:','HandleVisibility','off');
yline(ax,-hp.direction_confidence2_min,'k:','HandleVisibility','off');
for j=1:hp.N
    if isfinite(out.direction_lock_time_s(j))
        xline(ax,out.direction_lock_time_s(j),':','Color',S.USV(j,:), ...
            'LineWidth',0.8,'HandleVisibility','off');
    end
end
xlabel(ax,'physical time (s)');
ylabel(ax,'layer-2 evidence $\eta_2^{(j)}$');
title(ax,'(d) Finite direction commissioning','FontWeight','normal');
xlim(ax,[0 2]); local_grid(ax);
text(ax,0.97,0.08,sprintf('locked by %.1f s',max(out.direction_lock_time_s)), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',7.0,'BackgroundColor','w');

S.export(fig,fig_dir,'fig_hydro_closedloop');
fprintf('  fig_hydro_closedloop\n');
end

function local_grid(ax)
grid(ax,'on'); ax.GridLineStyle=':'; ax.GridAlpha=0.3;
end
