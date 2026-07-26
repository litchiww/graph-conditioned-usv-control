function usv_plots_app(results, mc_results, fig_dir, OPT)

if nargin < 4 || isempty(OPT)
    OPT.mc_mode = 'relabel';
    OPT.traj_rms = 'steady';
end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
S = usv_style(); S.apply();
cfg = usv_config(); layout = cfg.scenario_layout();

fig_A1(results, fig_dir, S, layout, OPT);  fprintf('  fig_trajectory\n');
fig_A2(results, fig_dir, S, layout);       fprintf('  fig_separation\n');
fig_A4(results, fig_dir, S, layout);       fprintf('  fig_seastate\n');
if ~strcmpi(OPT.mc_mode,'cut')
    fig_A5(mc_results, fig_dir, S, OPT); fprintf('  fig_montecarlo\n');
end
fig_A6(results, fig_dir, S);               fprintf('  fig_ablation_qbeta\n');
fig_A7(results, fig_dir, S);               fprintf('  fig_ablation_chi\n');
end

function fig_A1(R, fig_dir, S, layout, OPT)
cases = {'C0_baseline', '(a) C0: degree-weighted reference';
         'C2_full',     '(b) C2: graph-conditioned full method'};
tools = usv_metrics();

fig = S.figure('double',0.56);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

for k = 1:2
    a = R.(cases{k,1}).arr;
    ax = nexttile(tl); hold(ax,'on');

    for r = 1:size(layout.turbines,1)
        xc = layout.turbines(r,1); yc = layout.turbines(r,2);
        rectangle(ax,'Position',[xc-layout.safety_radius, yc-layout.safety_radius, ...
            2*layout.safety_radius, 2*layout.safety_radius], ...
            'Curvature',[1 1],'EdgeColor',[0.55 0.55 0.55],'LineStyle',':','LineWidth',0.55);
        plot(ax,xc,yc,'o','Color',[0.27 0.27 0.27],'MarkerFaceColor',[0.27 0.27 0.27],'MarkerSize',4);
        plot(ax,[xc xc],[yc-25 yc+25],'Color',[0.27 0.27 0.27],'LineWidth',0.9);
        plot(ax,[xc-35 xc+35],[yc yc],'Color',[0.4 0.4 0.4],'LineWidth',0.8);
    end
    for yc = layout.lane_centers_y
        plot(ax,layout.x_corridor,[yc yc],'--','Color',[0.6 0.6 0.6],'LineWidth',0.55);
    end

    lane_dev_max = 0;
    for j=1:3
        [X,Y,psi] = local_phys(a.t,a.x1(:,j),a.x2(:,j),j,layout);
        plot(ax,X,Y,S.USVLineStyle{j}, ...
            'Color',S.USV(j,:),'LineWidth',1.5);
        lane_dev_max = max(lane_dev_max, max(abs(Y - layout.lane_centers_y(j))));
        for ts = [0 10 30 60]
            [~,ii] = min(abs(a.t - ts));
            local_draw_triangle(ax,X(ii),Y(ii),psi(ii),55,S.USV(j,:));
        end
    end

    for ts = [0 10 30 60]
        [~,ii] = min(abs(a.t - ts));
        sync_x = zeros(3,1); sync_y = zeros(3,1);
        for jj = 1:3
            [Xj,Yj,~] = local_phys(a.t(ii),a.x1(ii,jj),a.x2(ii,jj),jj,layout);
            sync_x(jj) = Xj; sync_y(jj) = Yj;
        end
        plot(ax,sync_x,sync_y,'--','Color',[0.4 0.4 0.4],'LineWidth',0.8,'HandleVisibility','off');
        plot(ax,sync_x,sync_y,'.','Color',[0.3 0.3 0.3],'MarkerSize',6,'HandleVisibility','off');
    end

    [Xt,~,~] = local_phys(a.t,a.x1(:,3),a.x2(:,3),3,layout);
    for ts = [0 10 30 60]
        [~,ii] = min(abs(a.t - ts));
        text(ax,Xt(ii),1050,sprintf('t=%ds',ts),'HorizontalAlignment','center','FontSize',7.3);
    end

    tm = tools.task_metrics(a, layout);
    if strcmpi(OPT.traj_rms,'steady')
        rms_dev = tm.path_rms_m;
        window_txt = 'steady RMS';
    else
        rms_dev = sqrt(mean(a.x1(:).^2)) * layout.length_scale;
        window_txt = 'full RMS';
    end
    text(ax,0.97,0.92,sprintf('%s = %.1f m, max = %.0f m',window_txt,rms_dev,lane_dev_max), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',7.5,'Color',S.ink);

    xlabel(ax,'East (m)'); ylabel(ax,'North (m)');
    title(ax,cases{k,2},'FontWeight','normal');
    axis(ax,'equal');
    xlim(ax,layout.x_corridor); ylim(ax,layout.y_corridor);
    local_grid(ax);
end
S.export(fig, fig_dir, 'fig_trajectory');
end

function fig_A2(R, fig_dir, S, layout)
tools = usv_metrics();
a = R.C2_full.arr;
d = tools.inter_distance(a, layout);
fig = S.figure('double',0.47);
ax = axes('Parent',fig); hold(ax,'on');
patch(ax,[a.t(1) a.t(end) a.t(end) a.t(1)], ...
    [0 0 layout.usv_min_separation layout.usv_min_separation], ...
    [1.0 0.85 0.85],'EdgeColor','none','FaceAlpha',0.45, ...
    'DisplayName',sprintf('unsafe (<%g m)',layout.usv_min_separation));
plot(ax,d.t,d.d12,'-','Color',S.USV(1,:),'LineWidth',1.1,'DisplayName','USV1 <-> USV2');
plot(ax,d.t,d.d13,'--','Color',S.USV(2,:),'LineWidth',1.1,'DisplayName','USV1 <-> USV3');
plot(ax,d.t,d.d23,'-.','Color',S.USV(3,:),'LineWidth',1.1,'DisplayName','USV2 <-> USV3');
yline(ax,layout.usv_min_separation,'r:','LineWidth',0.8,'HandleVisibility','off');
dmin = min([d.d12; d.d13; d.d23]);
text(ax,0.985,0.09,sprintf('min separation = %.0f m (margin %.0f m)', ...
    dmin,dmin-layout.usv_min_separation), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
xlabel(ax,'time (s)'); ylabel(ax,'inter-USV distance (m)');
title(ax,'Inter-USV separation','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',4,'FontSize',6.8); local_grid(ax);
xlim(ax,[0 a.t(end)]);
S.export(fig, fig_dir, 'fig_separation');
end

function fig_A4(R, fig_dir, S, layout)
ss_keys = {'SS1','offshore_ss1';
           'SS2','offshore_nominal';
           'SS3','offshore_ss3';
           'SS4','offshore_rough'};
ss_labels = {'SS1 calm','SS2 nominal','SS3 moderate','SS4 rough'};
have = false(1,4); path_rms = nan(1,4); max_rat = nan(1,4); effort = nan(1,4);
for s = 1:4
    fld = sprintf('SS_%s', ss_keys{s,1});
    if isfield(R, fld)
        a = R.(fld).arr; o = R.(fld).out;
        mask = a.t >= a.t(end)/2;
        path_rms(s) = sqrt(mean(a.x1(mask,:).^2,'all')) * layout.length_scale;
        max_rat(s) = o.max_ratio;
        effort(s) = trapz(a.t, sum(a.u.^2,2));
        have(s) = true;
    end
end
sel = find(have);
if isempty(sel), fprintf('  (no SS_* results found; skipping A4)\n'); return; end
colors = [0.46 0.70 0.85; 0.37 0.68 0.37; 0.95 0.78 0.36; 0.90 0.50 0.30];

fig = S.figure('double',0.55);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); b = bar(ax,1:numel(sel),path_rms(sel),'FaceColor','flat','HandleVisibility','off');
b.CData = colors(sel,:);
set(ax,'XTick',1:numel(sel),'XTickLabel',ss_labels(sel));
xtickangle(ax,30);
ylabel(ax,'Lane-offset RMS (m)'); title(ax,'(a) Tracking','FontWeight','normal'); local_grid(ax); ax.YGrid='on';
for j=1:numel(sel), text(ax,j,path_rms(sel(j))*1.04,sprintf('%.3g',path_rms(sel(j))),'HorizontalAlignment','center','FontSize',7); end

ax = nexttile(tl); b = bar(ax,1:numel(sel),max_rat(sel),'FaceColor','flat','HandleVisibility','off');
b.CData = colors(sel,:);
yline(ax,1.0,'r--','LineWidth',0.8,'DisplayName','violation bound');
set(ax,'XTick',1:numel(sel),'XTickLabel',ss_labels(sel));
xtickangle(ax,30);
ylabel(ax,'max $|x_2|/F$'); title(ax,'(b) Envelope','FontWeight','normal');
legend(ax,'show','Location','northwest','FontSize',7); local_grid(ax); ax.YGrid='on';
ylim(ax,[0 max(1.15,max(max_rat(sel))*1.18)]);
for j=1:numel(sel), text(ax,j,max_rat(sel(j))+0.035,sprintf('%.3g',max_rat(sel(j))),'HorizontalAlignment','center','FontSize',7); end

ax = nexttile(tl); b = bar(ax,1:numel(sel),effort(sel),'FaceColor','flat','HandleVisibility','off');
b.CData = colors(sel,:);
set(ax,'XTick',1:numel(sel),'XTickLabel',ss_labels(sel));
xtickangle(ax,30);
ylabel(ax,'Normalized effort $\int\!\sum_j u_j^2dt$'); title(ax,'(c) Effort','FontWeight','normal'); local_grid(ax); ax.YGrid='on';
for j=1:numel(sel), text(ax,j,effort(sel(j))*1.04,sprintf('%.3g',effort(sel(j))),'HorizontalAlignment','center','FontSize',7); end

S.export(fig, fig_dir, 'fig_seastate');
end

function fig_A5(mc, fig_dir, S, ~)
if isempty(mc) || isempty(fieldnames(mc))
    fprintf('  (no MC results passed; skipping A5)\n'); return;
end
if ~isfield(mc,'group')
    for i=1:numel(mc), mc(i).group = 'combined_stress'; end
end

groups = {'in_domain','sign_only','gain_only','initial_only','load_only','combined_stress'};
labels = {'verified','signs','gains','initial','load','combined'};
keep = cellfun(@(g)any(strcmp({mc.group},g)),groups);
groups = groups(keep); labels = labels(keep);
outcomes = zeros(numel(groups),4);
completed = cell(numel(groups),1);
for gi=1:numel(groups)
    sub = mc(strcmp({mc.group},groups{gi}));
    [outcomes(gi,:), completed{gi}] = local_mc_summary(sub);
end

fig = S.figure('double',0.52);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
cols = [0.27 0.62 0.33; 0.90 0.65 0.41; 0.80 0.40 0.40; 0.45 0.45 0.45];
group_cols = lines(max(6,numel(groups)));

ax = nexttile(tl); hold(ax,'on');
bh = bar(ax,outcomes,'stacked','EdgeColor','none');
for k=1:4, bh(k).FaceColor=cols(k,:); end
set(ax,'XTick',1:numel(groups),'XTickLabel',labels); xtickangle(ax,25);
ylabel(ax,'runs'); title(ax,'(a) Factor-isolated outcomes','FontWeight','normal');
legend(ax,{'track + envelope','envelope, marginal track','envelope termination','numerical/exception'}, ...
    'Location','southoutside','NumColumns',2,'FontSize',6.2);
local_grid(ax); ax.YGrid='on';

ax = nexttile(tl); hold(ax,'on');
for gi=1:numel(groups)
    sub = mc(strcmp({mc.group},groups{gi}));
    ix = completed{gi};
    if isempty(ix), continue; end
    scatter(ax,[sub(ix).path_rms_m],[sub(ix).max_ratio],24,group_cols(gi,:), ...
        'filled','MarkerEdgeColor','k','DisplayName',labels{gi});
end
yline(ax,1,'r--','LineWidth',0.9,'HandleVisibility','off');
xlabel(ax,'path RMS error (m)'); ylabel(ax,'max $|x_2|/F$');
title(ax,'(b) Completed-run metrics','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.2);
local_grid(ax);

S.export(fig, fig_dir, 'fig_montecarlo');
end

function [parts, completed_ix] = local_mc_summary(mc)
if isempty(mc), parts=zeros(1,4); completed_ix=[]; return; end
reason = {mc.fail_reason};
env = strcmp(reason,'constraint_violation');
num = strcmp(reason,'numerical_limit') | strcmp(reason,'exception');
completed_ix = find(~env & ~num);
paper_ok = false(size(reason)); safe_ok = false(size(reason));
paper_ok(completed_ix) = [mc(completed_ix).paper_ok];
safe_ok(completed_ix) = [mc(completed_ix).safe_ok];
parts = [sum(paper_ok), sum(safe_ok & ~paper_ok), sum(env), sum(num)];
end

function fig_A6(R, fig_dir, S)
if ~isfield(R,'C8_raw_error'), fprintf('  (C8_raw_error not found; skipping A6)\n'); return; end
tools = usv_metrics();
a2 = R.C2_full.arr; o2 = R.C2_full.out;
a8 = R.C8_raw_error.arr; o8 = R.C8_raw_error.out;
fig = S.figure('double',0.60);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on');
semilogy(ax,a2.t,tools.moving_rms(vecnorm(a2.e,2,2),a2.t,5.0),'-','Color',S.safee,'LineWidth',1.4);
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'5-s moving RMS of $\|e\|$'); title(ax,'(a) C2 graph-conditioned coordinate','FontWeight','normal');
text(ax,0.96,0.90,sprintf('$\\mathrm{RMS}_{ss}=%.3f$',o2.rms_e_steady),'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
local_grid(ax); xlim(ax,[0 a2.t(end)]); ylim(ax,[1e-2 3]);

ax = nexttile(tl); hold(ax,'on');
semilogy(ax,a8.t,tools.moving_rms(vecnorm(a8.e,2,2),a8.t,5.0),'-','Color',[0.831 0.267 0.267],'LineWidth',1.4);
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'5-s moving RMS of $\|e\|$'); title(ax,'(b) C8 raw-error coordinate','FontWeight','normal');
text(ax,0.96,0.90,sprintf('$\\mathrm{RMS}_{ss}=%.3f$',o8.rms_e_steady),'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
local_grid(ax); xlim(ax,[0 a8.t(end)]); ylim(ax,[1e-2 3]);

ax = nexttile(tl); hold(ax,'on');
ratio2 = abs(a2.x2)./a2.F;
for j=1:3
    plot(ax,a2.t,ratio2(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.95,'DisplayName',sprintf('USV %d',j));
end
yline(ax,1.0,'r--','LineWidth',0.8,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'$|x_2^j|/F_2^j$'); title(ax,'(c) C2 full: envelope ratio','FontWeight','normal');
text(ax,0.04,0.88,sprintf('max=%.3f<1',o2.max_ratio), ...
    'Units','normalized','HorizontalAlignment','left','FontSize',7.5, ...
    'BackgroundColor','w','Margin',1.5);
legend(ax,'show','Location','northeast','FontSize',6.5); local_grid(ax); xlim(ax,[0 a2.t(end)]); ylim(ax,[0 1.05]);
local_ratio_inset(ax,a2.t,ratio2,S,[0 a2.t(end)],[0 0.08],[0.56 0.16 0.36 0.32]);

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a8.t,abs(a8.x2(:,j))./a8.F(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.95,'DisplayName',sprintf('USV %d',j));
end
yline(ax,1.0,'r--','LineWidth',0.8,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'$|x_2^j|/F_2^j$'); title(ax,'(d) C8 raw-error envelope ratio','FontWeight','normal');
text(ax,0.04,0.88,sprintf('max=%.3f<1',o8.max_ratio), ...
    'Units','normalized','HorizontalAlignment','left','FontSize',7.5, ...
    'BackgroundColor','w','Margin',1.5);
legend(ax,'show','Location','northeast','FontSize',6.5); local_grid(ax); xlim(ax,[0 a8.t(end)]); ylim(ax,[0 1.05]);

S.export(fig, fig_dir, 'fig_ablation_qbeta');
end

function fig_A7(R, fig_dir, S)
if ~isfield(R,'C7_no_chi'), fprintf('  (C7_no_chi not found; skipping A7)\n'); return; end
tools = usv_metrics();
a2 = R.C2_full.arr; o2 = R.C2_full.out;
a7 = R.C7_no_chi.arr; o7 = R.C7_no_chi.out;
fig = S.figure('double',0.60);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on');
semilogy(ax,a2.t,tools.moving_rms(vecnorm(a2.e,2,2),a2.t,5.0),'-','Color',S.safee,'LineWidth',1.4);
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'5-s moving RMS of $\|e\|$'); title(ax,'(a) C2 barrier coordinate','FontWeight','normal');
text(ax,0.96,0.90,sprintf('$\\mathrm{RMS}_{ss}=%.3f$',o2.rms_e_steady),'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
local_grid(ax); xlim(ax,[0 a2.t(end)]); ylim(ax,[1e-2 3]);

ax = nexttile(tl); hold(ax,'on');
semilogy(ax,a7.t,tools.moving_rms(vecnorm(a7.e,2,2),a7.t,5.0),'-','Color',[0.886 0.416 0.416],'LineWidth',1.4);
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'5-s moving RMS of $\|e\|$'); title(ax,'(b) C7 without $\chi$ coordinate','FontWeight','normal');
text(ax,0.96,0.90,sprintf('$\\mathrm{RMS}_{ss}=%.3f$',o7.rms_e_steady),'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
local_grid(ax); xlim(ax,[0 a7.t(end)]); ylim(ax,[1e-2 3]);

ax = nexttile(tl); hold(ax,'on');
ratio2 = abs(a2.x2)./a2.F;
for j=1:3
    plot(ax,a2.t,ratio2(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.95,'DisplayName',sprintf('USV %d',j));
end
yline(ax,1.0,'r--','LineWidth',0.8,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'$|x_2^j|/F_2^j$'); title(ax,'(c) C2 full: envelope ratio','FontWeight','normal');
text(ax,0.96,0.90,sprintf('max=%.3f<1',o2.max_ratio),'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
legend(ax,'show','Location','northeast','FontSize',6.5); local_grid(ax); xlim(ax,[0 a2.t(end)]); ylim(ax,[0 1.05]);
local_ratio_inset(ax,a2.t,ratio2,S,[0 a2.t(end)],[0 0.08],[0.56 0.16 0.36 0.32]);

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a7.t,abs(a7.x2(:,j))./a7.F(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.95,'DisplayName',sprintf('USV %d',j));
end
yline(ax,1.0,'r--','LineWidth',0.8,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'$|x_2^j|/F_2^j$'); title(ax,'(d) C7 without $\chi$: envelope ratio','FontWeight','normal');
text(ax,0.96,0.90,sprintf('max=%.3f<1',o7.max_ratio),'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
legend(ax,'show','Location','northeast','FontSize',6.5); local_grid(ax); xlim(ax,[0 a7.t(end)]); ylim(ax,[0 1.05]);

S.export(fig, fig_dir, 'fig_ablation_chi');
end

function [X,Y,psi] = local_phys(t,x1,x2,row,layout)
Y = layout.lane_centers_y(row) + x1*layout.length_scale;
X = layout.path_progress_rate * t;
psi = atan2(x2, 1.0);
end

function local_draw_triangle(ax, x, y, psi, sz, color)
ang0 = psi + [pi/2, pi/2 + 2*pi/3, pi/2 + 4*pi/3];
xs = x + sz*cos(ang0); ys = y + sz*sin(ang0);
patch(ax,xs,ys,color,'EdgeColor','k','LineWidth',0.6);
end

function local_grid(ax)
grid(ax,'on');
ax.GridLineStyle = ':';
ax.GridAlpha = 0.30;
end

function axI = local_inset_axes(ax, relpos)
fig = ancestor(ax,'figure');
p = ax.Position;
pos = [p(1)+relpos(1)*p(3), p(2)+relpos(2)*p(4), relpos(3)*p(3), relpos(4)*p(4)];
axI = axes('Parent',fig,'Position',pos);
hold(axI,'on');
set(axI,'Box','on','FontSize',6.5,'LineWidth',0.55);
end

function local_ratio_inset(ax,t,ratio,S,xr,yr,relpos)
axI = local_inset_axes(ax,relpos);
for j=1:size(ratio,2)
    plot(axI,t,ratio(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.55,'HandleVisibility','off');
end
xlim(axI,xr); ylim(axI,yr);
set(axI,'XTick',xr,'YTick',yr);
grid(axI,'on'); axI.GridLineStyle=':'; axI.GridAlpha=0.20;
end
