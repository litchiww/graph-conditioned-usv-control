function usv_plots_control(results, fig_dir, OPT)

if nargin < 3 || ~isfield(OPT,'keep_overview_consensus')
    OPT.keep_overview_consensus = true;
end

if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
S = usv_style(); S.apply();

fig_F1(results, fig_dir, S);  fprintf('  fig_main\n');
fig_F4(results, fig_dir, S);  fprintf('  fig_topology\n');
fig_F5(results, fig_dir, S);  fprintf('  fig_decaying_diagnostic\n');
if OPT.keep_overview_consensus
    fig_F8(results, fig_dir, S);  fprintf('  fig_consensus\n');
    fig_F9(results, fig_dir, S);  fprintf('  fig_overview\n');
end
fig_search_adaptive(results, fig_dir, S);  fprintf('  fig_search_adaptive\n');
fig_nussbaum(results, fig_dir, S);         fprintf('  fig_nussbaum\n');
end

function fig_F1(R, fig_dir, S)

a = R.C2_full.arr;
o = R.C2_full.out;
tools = usv_metrics();
fig = S.figure('double', 0.68);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on');
patch(ax,[a.t(end)/2 a.t(end) a.t(end) a.t(end)/2], ...
      [-10 -10 10 10],[0.85 0.95 0.85], ...
      'EdgeColor','none','FaceAlpha',0.18,'HandleVisibility','off');
for j=1:3
    plot(ax,a.t,a.x1(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.3, ...
        'DisplayName',sprintf('USV %d',j));
end
plot(ax,a.t,a.yd,'k--','LineWidth',1.0,'DisplayName','$y_d$');
xlabel(ax,'time (s)'); ylabel(ax,'$x_1$ and $y_d$'); title(ax,'(a) Cooperative path-following','FontWeight','normal');
ylim(ax,[-0.8 0.5]); xlim(ax,[0 a.t(end)]);
legend(ax,'show','Location','northoutside','NumColumns',4,'FontSize',7); local_grid(ax);

ax = nexttile(tl); hold(ax,'on');
en = tools.moving_rms(vecnorm(a.e,2,2), a.t, 1.0);
qn = tools.moving_rms(vecnorm(a.q,2,2), a.t, 1.0);
cn = tools.moving_rms(vecnorm(a.chi,2,2), a.t, 1.0);
semilogy(ax,a.t,en,'-','LineWidth',1.2,'DisplayName','$\|e\|$');
semilogy(ax,a.t,qn,'--','LineWidth',1.2,'DisplayName','$\|q_\beta\|$');
semilogy(ax,a.t,cn,'-.','LineWidth',1.1,'Color',[0.85 0.65 0.13],'DisplayName','$\|\chi\|$');
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'norm (1-s moving RMS, log)');
title(ax,'(b) Closed-loop norm histories','FontWeight','normal');
legend(ax,'show','Location','northeast','FontSize',7); local_grid(ax);
xlim(ax,[0 a.t(end)]); ylim(ax,[8e-4 5]);

ax = nexttile(tl); hold(ax,'on');
ratio = abs(a.x2)./a.F;
for j=1:3
    plot(ax,a.t,ratio(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.0, ...
        'DisplayName',sprintf('USV %d',j));
end
yline(ax,1.0,'r--','LineWidth',0.8,'DisplayName','bound');
xlabel(ax,'time (s)'); ylabel(ax,'$|x_2^j|/F_2^j$'); title(ax,'(c) Operational-envelope ratio','FontWeight','normal');
text(ax,0.97,0.90,sprintf('max = %.3f; bound = 1',o.max_ratio), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.6); local_grid(ax);
xlim(ax,[0 a.t(end)]); ylim(ax,[0 1.05]);
local_ratio_inset(ax,a.t,ratio,S,[0 a.t(end)],[0 0.08],[0.56 0.16 0.36 0.32]);

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a.t,a.u(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.95, ...
        'DisplayName',sprintf('$u_%d$',j));
end
xlabel(ax,'time (s)'); ylabel(ax,'control command'); title(ax,'(d) Control inputs','FontWeight','normal');
text(ax,0.97,0.90,sprintf('max |u| = %.2f',o.max_abs_u), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7.5);
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.6); local_grid(ax);
xlim(ax,[0 a.t(end)]);
m = 1.05*max(abs(a.u(:)));
if m == 0, m = 1; end
ylim(ax,[-m m]);
local_signal_inset(ax,a.t,a.u,S,[8 a.t(end)],[-0.5 0.5],[0.10 0.16 0.38 0.32]);

S.export(fig, fig_dir, 'fig_main');
end

function fig_F4(R, fig_dir, S)
tools = usv_metrics();
topo = {'C2_full','chain','-',  [0.122 0.467 0.706];
        'C3_star','star','--', [1.000 0.498 0.055];
        'C4_complete','complete','-.', [0.173 0.627 0.173];
        'C5_directed','directed',':', [0.839 0.153 0.157]};
fig = S.figure('double',0.52);
ax = axes('Parent',fig); hold(ax,'on');
T_end = R.C2_full.arr.t(end);
patch(ax,[T_end/2 T_end T_end T_end/2],[1e-3 1e-3 1e3 1e3], ...
    [0.85 0.95 0.85],'EdgeColor','none','FaceAlpha',0.18,'DisplayName','steady window');
for k=1:size(topo,1)
    a = R.(topo{k,1}).arr;
    rms = tools.moving_rms(vecnorm(a.e,2,2),a.t,5.0);
    semilogy(ax,a.t,rms,topo{k,3},'Color',topo{k,4},'LineWidth',1.2,'DisplayName',topo{k,2});
end
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'5-s moving RMS of $\|e\|$');
title(ax,'Topology comparison','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',5,'FontSize',7.2); local_grid(ax);
xlim(ax,[0 T_end]); ylim(ax,[1e-2 2]);
S.export(fig, fig_dir, 'fig_topology');
end

function fig_F5(R, fig_dir, S)
a = R.C9_decaying_diagnostic.arr;
tools = usv_metrics();
fig = S.figure('double',0.48);
ax = axes('Parent',fig); hold(ax,'on');
flo = 1e-6;
en = max(tools.moving_rms(vecnorm(a.e,2,2),   a.t, 1.0), flo);
qn = max(tools.moving_rms(vecnorm(a.q,2,2),   a.t, 1.0), flo);
cn = max(tools.moving_rms(vecnorm(a.chi,2,2), a.t, 1.0), flo);
semilogy(ax,a.t,en,'-','DisplayName','$\|e\|$');
semilogy(ax,a.t,qn,'--','DisplayName','$\|q_\beta\|$');
semilogy(ax,a.t,cn,'-.','Color',[0.85 0.65 0.13],'DisplayName','$\|\chi\|$');
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'norm (1-s moving RMS, log)');
title(ax,'Decaying-disturbance diagnostic','FontWeight','normal');
text(ax,0.04,0.10,'late regrowth prevents an asymptotic claim', ...
    'Units','normalized','HorizontalAlignment','left','FontSize',7.5);
legend(ax,'show','Location','northeast','FontSize',7.5); local_grid(ax);
xlim(ax,[0 150]); ylim(ax,[flo 2]);
S.export(fig, fig_dir, 'fig_decaying_diagnostic');
end

function fig_F8(R, fig_dir, S)
a = R.C2_full.arr;
tools = usv_metrics();
fig = S.figure('double',0.50);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on');
diffs = [a.x1(:,1)-a.x1(:,2), a.x1(:,2)-a.x1(:,3), a.x1(:,1)-a.x1(:,3)];
plot(ax,a.t,diffs(:,1),'-','Color',S.USV(1,:),'LineWidth',1.1,'DisplayName','$x_1^1-x_1^2$');
plot(ax,a.t,diffs(:,2),'--','Color',S.USV(2,:),'LineWidth',1.1,'DisplayName','$x_1^2-x_1^3$');
plot(ax,a.t,diffs(:,3),'-.','Color',S.USV(3,:),'LineWidth',1.1,'DisplayName','$x_1^1-x_1^3$');
yline(ax,0,'k:','LineWidth',0.6,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'pairwise difference'); title(ax,'(a) Pairwise consensus differences','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.5); local_grid(ax); xlim(ax,[0 a.t(end)]);
md = 1.05*max(abs(diffs(:)));
if md == 0, md = 1; end
ylim(ax,[-md md]);

ax = nexttile(tl); hold(ax,'on');
d12 = abs(a.x1(:,1)-a.x1(:,2));
d23 = abs(a.x1(:,2)-a.x1(:,3));
d13 = abs(a.x1(:,1)-a.x1(:,3));
semilogy(ax,a.t,tools.moving_rms(d12,a.t,3.0),'-','Color',S.USV(1,:),'LineWidth',1.1,'DisplayName','$|x_1^1-x_1^2|$ RMS');
semilogy(ax,a.t,tools.moving_rms(d23,a.t,3.0),'--','Color',S.USV(2,:),'LineWidth',1.1,'DisplayName','$|x_1^2-x_1^3|$ RMS');
semilogy(ax,a.t,tools.moving_rms(d13,a.t,3.0),'-.','Color',S.USV(3,:),'LineWidth',1.1,'DisplayName','$|x_1^1-x_1^3|$ RMS');
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'pairwise |difference| (3-s RMS, log)');
title(ax,'(b) Convergence rate','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.5); local_grid(ax); xlim(ax,[0 a.t(end)]);
S.export(fig, fig_dir, 'fig_consensus');
end

function fig_F9(R, fig_dir, S)
tools = usv_metrics();
cs = S.ctrl_styles();
fig = S.figure('double',0.48);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
T_end = R.C2_full.arr.t(end);

ax = nexttile(tl); hold(ax,'on');
for k=1:size(cs,1)
    a = R.(cs{k,1}).arr;
    r = tools.moving_rms(vecnorm(a.e,2,2),a.t,5.0);
    semilogy(ax,a.t,r,cs{k,3},'Color',cs{k,5},'LineWidth',cs{k,4},'DisplayName',cs{k,2});
end
set(ax,'YScale','log');
xlabel(ax,'time (s)'); ylabel(ax,'5-s moving RMS of $\|e\|$'); title(ax,'(a) Tracking','FontWeight','normal');
lgd = legend(ax,'show','Orientation','horizontal','NumColumns',5,'FontSize',6.3);
lgd.Layout.Tile = 'south';
local_grid(ax); ylim(ax,[1e-2 2]); xlim(ax,[0 T_end]);

ax = nexttile(tl); hold(ax,'on');
for k=1:size(cs,1)
    a = R.(cs{k,1}).arr;
    ratio = max(abs(a.x2)./a.F,[],2);
    plot(ax,a.t,ratio,cs{k,3},'Color',cs{k,5},'LineWidth',cs{k,4},'DisplayName',cs{k,2});
end
yline(ax,1.0,'r--','LineWidth',0.8,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'$\max_j |x_2^j|/F_2^j$'); title(ax,'(b) Envelope ratio','FontWeight','normal');
local_grid(ax); ylim(ax,[0 1.05]); xlim(ax,[0 T_end]);

ax = nexttile(tl); hold(ax,'on');
pk = zeros(size(cs,1),1);
for k=1:size(cs,1)
    pk(k) = max(abs(R.(cs{k,1}).arr.u(:)));
end
b = bar(ax,1:size(cs,1),pk,'FaceColor','flat','EdgeColor',[0.2 0.2 0.2], ...
    'LineWidth',0.5,'HandleVisibility','off');
for k=1:size(cs,1), b.CData(k,:) = cs{k,5}; end
set(ax,'YScale','log','XTick',1:size(cs,1), ...
    'XTickLabel',{'C0','C1','C2','C10','C11'},'TickLabelInterpreter','latex');
ylabel(ax,'peak $|u|$ (log)'); title(ax,'(c) Control effort (peak)','FontWeight','normal');
local_grid(ax); ax.YGrid = 'on'; ylim(ax,[1 max(pk)*1.8]);
for k=1:size(cs,1)
    text(ax,k,pk(k)*1.12,sprintf('%.1f',pk(k)), ...
        'HorizontalAlignment','center','FontSize',7.5,'FontWeight','bold');
end
S.export(fig, fig_dir, 'fig_overview');
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

function local_signal_inset(ax,t,y,S,xr,yr,relpos)
axI = local_inset_axes(ax,relpos);
for j=1:size(y,2)
    plot(axI,t,y(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',0.55,'HandleVisibility','off');
end
yline(axI,0,'k:','LineWidth',0.45,'HandleVisibility','off');
xlim(axI,xr); ylim(axI,yr);
set(axI,'XTick',xr,'YTick',yr);
grid(axI,'on'); axI.GridLineStyle=':'; axI.GridAlpha=0.20;
end

function fig_search_adaptive(R, fig_dir, S)

a = R.C2_full.arr;
fig = S.figure('double', 0.68);
tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a.t,a.s1(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.25, ...
        'DisplayName',sprintf('$s_1^{%d}$',j));
end
xlabel(ax,'time (s)'); ylabel(ax,'$s_1^j$');
title(ax,'(a) Layer-1 search variable','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.5); local_grid(ax); xlim(ax,[0 a.t(end)]);

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a.t,a.s2(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.25, ...
        'DisplayName',sprintf('$s_2^{%d}$',j));
end
xlabel(ax,'time (s)'); ylabel(ax,'$s_2^j$');
title(ax,'(b) Layer-2 search variable','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.5); local_grid(ax); xlim(ax,[0 a.t(end)]);

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a.t,a.theta1(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.25, ...
        'DisplayName',sprintf('$\\hat\\theta_1^{%d}$',j));
end
xlabel(ax,'time (s)'); ylabel(ax,'$\hat\theta_1^j$ (estimate)');
title(ax,'(c) Layer-1 adaptive estimate','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.5); local_grid(ax); xlim(ax,[0 a.t(end)]);

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a.t,a.theta2(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.25, ...
        'DisplayName',sprintf('$\\hat\\theta_2^{%d}$',j));
end
xlabel(ax,'time (s)'); ylabel(ax,'$\hat\theta_2^j$ (estimate)');
title(ax,'(d) Layer-2 adaptive estimate','FontWeight','normal');
legend(ax,'show','Location','southoutside','NumColumns',3,'FontSize',6.5); local_grid(ax); xlim(ax,[0 a.t(end)]);
ylim(ax,[0 0.08]);

S.export(fig, fig_dir, 'fig_search_adaptive');
end

function fig_nussbaum(R, fig_dir, S)

a = R.C2_full.arr; p = R.C2_full.p;

g1 = [-1.2;  1.2; -0.8];
g2 = [ 2.0; -2.0;  1.5];

if isfield(a,'N1_effective') && isfield(a,'N2_effective')
    M1=a.N1_effective; M2=a.N2_effective;
else
    M1=local_nussbaum_eval(a.s1,p.eps1,p.nussbaum_index(1,:),p);
    M2=local_nussbaum_eval(a.s2,p.eps2,p.nussbaum_index(2,:),p);
end
gM1=M1.*g1.';
gM2=M2.*g2.';
tx=min(0.6,a.t(end));

fig = S.figure('double', 0.46);
tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on');
for j=1:3
    plot(ax,a.t,a.direction_blend2(:,j),S.USVLineStyle{j}, ...
        'Color',S.USV(j,:),'LineWidth',1.25, ...
        'DisplayName',sprintf('$w_2^{%d}$',j));
end
yline(ax,1,'k:','LineWidth',0.6,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'confidence blend $w_2^j$');
title(ax,'(a) Layer-2 local confidence blend','FontWeight','normal');
legend(ax,'show','Location','southeast','FontSize',7); local_grid(ax);
xlim(ax,[0 tx]); ylim(ax,[0 1.05]);

ax = nexttile(tl); hold(ax,'on');
mt=a.t<=tx;
allv = [gM1(mt,:); gM2(mt,:)]; allv=allv(:);
yl = [min(allv)*1.18, max(allv)*1.18];
patch(ax,[0 tx tx 0],[0 0 yl(1) yl(1)], ...
    [0.91 0.95 0.92],'EdgeColor','none','FaceAlpha',0.7, ...
    'DisplayName','stabilizing region ($gM<0$)');
for j=1:3
    plot(ax,a.t,gM2(:,j),'-','Color',S.USV(j,:),'LineWidth',1.3, ...
        'DisplayName',sprintf('$g_2^{%d}M_2^{%d}$',j,j));
    plot(ax,a.t,gM1(:,j),'--','Color',S.USV(j,:),'LineWidth',1.0, ...
        'HandleVisibility','off');
end
yline(ax,0,'k-','LineWidth',0.6,'HandleVisibility','off');
xlabel(ax,'time (s)'); ylabel(ax,'in-loop multiplier $g_i^jM_i^j$');
title(ax,'(b) In-loop multipliers enter the stabilizing region','FontWeight','normal');
text(ax,0.045,0.92,'unknown true signs differ', ...
    'Units','normalized','FontSize',6.6,'Color',[0.2 0.2 0.2], ...
    'BackgroundColor','w','Margin',1.5);
text(ax,0.045,0.84,'$g_1=[-,+,-],\; g_2=[+,-,+]$', ...
    'Units','normalized','FontSize',6.6,'Color',[0.2 0.2 0.2], ...
    'BackgroundColor','w','Margin',1.5);
text(ax,0.045,0.76,'$M_i=(1-w_i)N_i+w_i\hat d_i$', ...
    'Units','normalized','FontSize',6.6,'Color',[0.2 0.2 0.2], ...
    'BackgroundColor','w','Margin',1.5);
legend(ax,'show','Location','southoutside','NumColumns',2,'FontSize',6.4); local_grid(ax);
xlim(ax,[0 tx]); ylim(ax,yl);

S.export(fig, fig_dir, 'fig_nussbaum');
end

function N = local_nussbaum_eval(s, eps_, idx, p)

if strcmpi(p.nussbaum_family,'huang2018')
    alpha=p.nussbaum_alpha; beta=p.nussbaum_beta;
    w=beta.^(-idx); C=sqrt(alpha^2*beta.^(2*idx)+1)./beta.^idx;
    N=eps_*exp(min(alpha*abs(s),p.exp_clip)).*sin(s.*w).*C;
else
    k = p.nussbaum_base.^(idx-1);
    s2 = min(s.^2, p.exp_clip);
    N = eps_ * (2.*s.*sinh(s2).*cos(k.*s) ...
        - k.*(cosh(s2)-1).*sin(k.*s));
end
end
