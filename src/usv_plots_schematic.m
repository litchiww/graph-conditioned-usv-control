function usv_plots_schematic(~, case_dir, fig_dir)

if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
fig_architecture(fig_dir);
fig_concept(case_dir, fig_dir);
end

function fig_architecture(fig_dir)
S = usv_style(); S.apply();
fig = S.figure('double',0.55);
ax = axes('Parent',fig,'Position',[0 0 1 1]);
axis(ax,[0 100 0 70]); axis(ax,'off'); hold(ax,'on');

box = @(cx,cy,w,h,txt,fc,ec,fs) draw_box(ax,cx,cy,w,h,txt,fc,ec,fs);
arr = @(p0,p1) quiverarrow(ax,p0,p1,S.ink,1.0);

box(31,65,46,6.5,'Virtual leader: bounded reference $y_d(t)$',[1 0.965 0.878],[0.788 0.635 0.153],8);
draw_frame(ax,5,5.5,52,56,[0.541 0.592 0.651]);
text(ax,7.5,59,'Controller (vehicle j)','FontSize',7.5,'Color',[0.204 0.275 0.353]);

box(31,53,46,6.2,'Leader-following error  $e=Hx_1-by_d$',S.blockf,S.blocke,7.1);
box(31,44,46,7.0,{'Graph-conditioned coordinate','$q_\beta=D_\beta^{-1}H^\top e$'},S.hifill,S.hilite,8);
box(31,34.5,46,7.6,{'Layer 1 virtual control  $\alpha_1$','local Nussbaum search + adaptive estimate'},S.blockf,S.blocke,7.2);
box(31,24,46,7.8,{'Barrier coordinate / envelope','$\chi_2$ protects $|x_2|<F_2(t)$'},S.hifill,S.hilite,7.2);
box(31,14,46,7.4,{'Layer 2 real control  u','gradient self-cancel + local Nussbaum search'},S.blockf,S.blocke,7.2);
box(31,7.7,46,4.4,'Constant gains + smooth nonlinear damping', ...
    [0.941 0.925 0.969],[0.478 0.357 0.651],6.9);

arr([31 61.7],[31 56.2]); arr([31 49.9],[31 47.6]);
arr([31 40.4],[31 38.4]); arr([31 30.6],[31 27.9]);
arr([31 20.1],[31 17.8]); arr([31 10.3],[31 9.7]);
text(ax,55,44,'(I)','Color',S.hilite,'FontWeight','bold','FontSize',8);
text(ax,55,24,'(II)','Color',S.hilite,'FontWeight','bold','FontSize',8);

draw_frame(ax,64,14,32,40,[0.243 0.435 0.627]);
text(ax,80,50.5,'USV fleet plant  j = 1,...,N', ...
    'FontWeight','bold','FontSize',8,'Color',[0.153 0.290 0.420], ...
    'HorizontalAlignment','center');
text(ax,80,47.2,'network-coordinate idea','FontSize',6.9, ...
    'Color',[0.153 0.290 0.420],'HorizontalAlignment','center');
draw_mini_network(ax,S,71.5,43.5,'coupled');
draw_mini_arrow(ax,[77.0 43.5],[83.0 43.5],S.ink,0.75);
draw_mini_network(ax,S,88.5,43.5,'matched');
text(ax,71.5,39.2,{'entangled','searches'},'FontSize',5.8,'HorizontalAlignment','center','Color',S.hilite);
text(ax,88.5,39.2,{'matched','local laws'},'FontSize',5.8,'HorizontalAlignment','center','Color',S.safee);
text(ax,80,32.0,{'state layer', ...
    '$x_1$ dynamics:  $g_1^j x_2^j+f_1^j$', ...
    '$x_2$ dynamics:  $g_2^j u^j+f_2^j$', ...
    '$g_i^j$ fixed, unknown, nonidentical'}, ...
    'HorizontalAlignment','center','FontSize',7.1);

arr([54 14],[64 25]); text(ax,60,21,'$u^j$','FontSize',7);
box(80,8.5,32,6.2,{'Lyapunov guarantee','bounded signals, invariant envelope, $e$ UUB'},S.safef,S.safee,7);
arr([80 14],[80 11.7]);
draw_line(ax,[96 98],[28 28],S.blocke,1.0);
draw_line(ax,[98 98],[28 3],S.blocke,1.0);
draw_line(ax,[98 8],[3 3],S.blocke,1.0);
draw_line(ax,[8 8],[3 53],S.blocke,1.0);
arr([8 50],[8.5 53]);
text(ax,44,1.6,'feedback states $x_1,x_2$','FontSize',7,'HorizontalAlignment','center','Color',S.blocke);

S.export(fig, fig_dir, 'fig_architecture');
fprintf('  schematic: fig_architecture\n');
end

function fig_concept(case_dir, fig_dir)
S = usv_style(); S.apply();
Sdat = load(fullfile(case_dir,'C2_full.mat')); a = Sdat.arr;
j = 2;
fig = S.figure('single',0.78);
ax = axes('Parent',fig); hold(ax,'on');
m = a.t <= 1.5;
fill(ax,[a.t(m); flipud(a.t(m))],[a.F(m,j); -flipud(a.F(m,j))],S.safef, ...
    'EdgeColor','none','DisplayName','$\pm F_2^{(2)}$ envelope');
plot(ax,a.t(m), a.F(m,j),'-','Color',S.safee,'LineWidth',1.0,'HandleVisibility','off');
plot(ax,a.t(m),-a.F(m,j),'-','Color',S.safee,'LineWidth',1.0,'HandleVisibility','off');
plot(ax,a.t(m), a.alpha1(m,j),':','Color',[0.45 0.45 0.45],'LineWidth',1.0,'DisplayName','$\alpha_1^{(2)}(t)$');
plot(ax,a.t(m), a.x2(m,j),'-','Color',S.USV(j,:),'LineWidth',1.25,'DisplayName','$x_2^{(2)}(t)$');
yyaxis(ax,'right');
plot(ax,a.t(m),a.s1(m,j),'--','Color',S.hilite,'LineWidth',1.05,'DisplayName','$s_1^{(2)}(t)$');
ylabel(ax,'search variable $s_1^{(2)}$'); ax.YColor = S.hilite;
yyaxis(ax,'left');
xlabel(ax,'time (s)'); ylabel(ax,'$\pm F_2^{(2)}$, $x_2^{(2)}$, $\alpha_1^{(2)}$');
title(ax,'Envelope response during local direction handoff','FontWeight','normal');
xlim(ax,[0 1.5]); grid(ax,'on'); ax.GridLineStyle=':'; ax.GridAlpha=0.3;
ymax=max([a.F(m,j); abs(a.alpha1(m,j)); abs(a.x2(m,j))]);
ylim(ax,[-1.08*ymax 1.08*ymax]);
legend(ax,'Location','southoutside','NumColumns',4,'FontSize',6.4);
S.export(fig, fig_dir, 'fig_concept');
fprintf('  schematic: fig_concept\n');
end

function draw_mini_network(ax,S,cx,cy,mode)
xy = [cx-3.4 cy+0.8; cx cy+0.8; cx+3.4 cy+0.8];
for i=1:2
    plot(ax,xy(i:i+1,1),xy(i:i+1,2),'-','Color',[0.45 0.45 0.45], ...
        'LineWidth',0.65,'HandleVisibility','off');
end
for j=1:3
    patch(ax,xy(j,1)+0.55*cos(linspace(0,2*pi,30)), ...
        xy(j,2)+0.55*sin(linspace(0,2*pi,30)), S.USV(j,:), ...
        'EdgeColor','k','LineWidth',0.4,'HandleVisibility','off');
end
switch mode
    case 'coupled'
        plot(ax,[xy(1,1) xy(2,1) xy(3,1) xy(1,1)], ...
            [cy-1.5 xy(2,2) cy-1.5 cy-1.5], '--', ...
            'Color',S.hilite,'LineWidth',0.75,'HandleVisibility','off');
    case 'matched'
        for j=1:3
            draw_mini_arrow(ax,[xy(j,1) cy+0.25],[xy(j,1) cy-1.4],S.safee,0.65);
        end
end
end

function draw_mini_arrow(ax,p0,p1,col,lw)
plot(ax,[p0(1) p1(1)],[p0(2) p1(2)],'-','Color',col,'LineWidth',lw,'HandleVisibility','off');
ang = atan2(p1(2)-p0(2),p1(1)-p0(1));
L = 0.45;
hx = p1(1) - L*cos(ang-pi/7); hy = p1(2) - L*sin(ang-pi/7);
gx = p1(1) - L*cos(ang+pi/7); gy = p1(2) - L*sin(ang+pi/7);
patch(ax,[p1(1) hx gx],[p1(2) hy gy],col,'EdgeColor',col,'HandleVisibility','off');
end

function draw_chain_nodes(ax,S)
xy = [0.30 0.62; 0.50 0.62; 0.70 0.62];
for i=1:2
    plot(ax,xy(i:i+1,1),xy(i:i+1,2),'-','Color',[0.45 0.45 0.45],'LineWidth',0.9);
end
for j=1:3
    patch(ax,xy(j,1)+0.045*cos(linspace(0,2*pi,50)), ...
        xy(j,2)+0.045*sin(linspace(0,2*pi,50)), S.USV(j,:), ...
        'EdgeColor','k','LineWidth',0.5);
    text(ax,xy(j,1),xy(j,2)-0.12,sprintf('N_1^%d(s_1^%d)',j,j), ...
        'HorizontalAlignment','center','FontSize',6.8);
end
end

function draw_box(ax,cx,cy,w,h,txt,fc,ec,fs)
rectangle(ax,'Position',[cx-w/2 cy-h/2 w h],'Curvature',[0.10 0.20], ...
    'FaceColor',fc,'EdgeColor',ec,'LineWidth',1.0);
text(ax,cx,cy,txt,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fs);
end

function draw_frame(ax,x,y,w,h,ec)
rectangle(ax,'Position',[x y w h],'Curvature',[0.04 0.08],'EdgeColor',ec,'LineWidth',1.1);
end

function draw_line(ax,xs,ys,col,lw)
plot(ax,xs,ys,'-','Color',col,'LineWidth',lw,'HandleVisibility','off');
end

function quiverarrow(ax,p0,p1,col,lw)
plot(ax,[p0(1) p1(1)],[p0(2) p1(2)],'-','Color',col,'LineWidth',lw,'HandleVisibility','off');
ang = atan2(p1(2)-p0(2),p1(1)-p0(1)); L = 1.3;
if max(abs(p1-p0)) <= 1, L = 0.035; end
hx = p1(1) - L*cos(ang-pi/7); hy = p1(2) - L*sin(ang-pi/7);
gx = p1(1) - L*cos(ang+pi/7); gy = p1(2) - L*sin(ang+pi/7);
patch(ax,[p1(1) hx gx],[p1(2) hy gy],col,'EdgeColor',col,'HandleVisibility','off');
end
