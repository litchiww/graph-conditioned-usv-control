function usv_plots_redesign(root_dir,fig_dir)

if nargin<1, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
if nargin<2, fig_dir=fullfile(root_dir,'figures'); end
data_dir=fullfile(root_dir,'data');
S=usv_style(); S.apply();

P=readtable(fullfile(data_dir,'sign_redesign_patterns.csv'));
M=readtable(fullfile(data_dir,'sign_redesign_mc.csv'));
I=readtable(fullfile(data_dir,'sign_direction_identification_audit.csv'));

fig=S.figure('double',0.70);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); hold(ax,'on');
plot(ax,P.pattern_id,P.rms_e,'-o','Color',S.USV(1,:), ...
    'MarkerFaceColor','w','MarkerSize',2.8,'LineWidth',0.9);
yline(ax,0.10,'--','Color',S.hilite,'LineWidth',0.9);
xlabel(ax,'sign-pattern ID'); ylabel(ax,'steady RMS $\|e\|$');
title(ax,'(a) Exhaustive fixed signs','FontWeight','normal');
xlim(ax,[0 63]); ylim(ax,[0 0.105]); local_grid(ax);

ax=nexttile(tl); hold(ax,'on');
plot(ax,P.pattern_id,P.max_ratio,'-o','Color',S.safee, ...
    'MarkerFaceColor','w','MarkerSize',2.8,'LineWidth',0.9);
yline(ax,0.98,'--','Color',S.hilite,'LineWidth',0.9);
xlabel(ax,'sign-pattern ID'); ylabel(ax,'max $|x_2^j|/F_2^j$');
title(ax,'(b) Continuous envelope audit','FontWeight','normal');
xlim(ax,[0 63]); ylim(ax,[0.75 1.0]); local_grid(ax);

ax=nexttile(tl); hold(ax,'on');
yyaxis(ax,'left');
plot(ax,M.run_id,M.rms_e,'o','Color',S.USV(1,:), ...
    'MarkerFaceColor',S.USV(1,:),'MarkerSize',3.1);
ylabel(ax,'steady RMS $\|e\|$'); ylim(ax,[0 0.10]);
yyaxis(ax,'right');
plot(ax,M.run_id,M.max_ratio,'s','Color',S.hilite, ...
    'MarkerFaceColor','w','MarkerSize',3.1);
ylabel(ax,'max ratio'); ylim(ax,[0.70 1.0]);
xlabel(ax,'combined-stress run');
title(ax,'(c) Random signs and perturbations','FontWeight','normal');
xlim(ax,[0.5 height(M)+0.5]); local_grid(ax);

ax=nexttile(tl); hold(ax,'on');
plot(ax,I.pattern_id,I.latest_activation_s,'-','Color',S.blocke, ...
    'LineWidth',1.0,'DisplayName','first activation');
plot(ax,I.pattern_id,I.latest_99pct_blend_s,'--','Color',S.hilite, ...
    'LineWidth',1.0,'DisplayName','99\% handoff');
xlabel(ax,'sign-pattern ID'); ylabel(ax,'latest channel time (s)');
title(ax,'(d) Local direction handoff','FontWeight','normal');
xlim(ax,[0 63]); ylim(ax,[0 0.40]); local_grid(ax);
legend(ax,'Location','southoutside','NumColumns',2,'FontSize',7.0);

S.export(fig,fig_dir,'fig_sign_redesign');
fprintf('  fig_sign_redesign\n');
end

function local_grid(ax)
grid(ax,'on'); ax.GridLineStyle=':'; ax.GridAlpha=0.30;
end
