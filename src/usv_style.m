function S = usv_style()

S.apply        = @apply_defaults;
S.figure       = @make_figure;
S.export       = @export_fig;
S.USV          = [0.000 0.447 0.698;
                  0.835 0.369 0.000;
                  0.000 0.620 0.451];
S.USVLineStyle = {'-','--','-.'};
S.ctrl_styles  = @controller_styles;
S.W1 = 3.35;
S.W2 = 7.16;
S.ink     = [0.10 0.10 0.10];
S.hilite  = [0.75 0.23 0.17];
S.hifill  = [0.984 0.925 0.918];
S.blockf  = [0.933 0.949 0.965];
S.blocke  = [0.290 0.353 0.416];
S.safef   = [0.910 0.953 0.918];
S.safee   = [0.173 0.627 0.173];
end

function apply_defaults()
set(groot, 'defaultAxesFontName','Times New Roman', ...
           'defaultTextFontName','Times New Roman', ...
           'defaultAxesFontSize',8.5, ...
           'defaultTextFontSize',8.5, ...
           'defaultAxesTitleFontSizeMultiplier',1.0, ...
           'defaultAxesLabelFontSizeMultiplier',1.0, ...
           'defaultAxesLineWidth',0.7, ...
           'defaultLineLineWidth',1.2, ...
           'defaultAxesTickLabelInterpreter','latex', ...
           'defaultTextInterpreter','latex', ...
           'defaultLegendInterpreter','latex', ...
           'defaultColorbarTickLabelInterpreter','latex', ...
           'defaultAxesGridLineStyle',':', ...
           'defaultAxesGridAlpha',0.30);
end

function fig = make_figure(width_key, aspect)
if nargin < 2 || isempty(aspect), aspect = 0.62; end
S = usv_style();
switch lower(width_key)
    case 'single'
        w = S.W1;
    case 'double'
        w = S.W2;
    otherwise
        w = S.W2;
end
fig = figure('Visible','off','Units','inches', ...
             'Position',[1 1 w w*aspect], 'Color','w');
set(fig,'PaperUnits','inches','PaperPositionMode','auto');
end

function export_fig(fig, fig_dir, name)
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
exportgraphics(fig, fullfile(fig_dir,[name '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig, fullfile(fig_dir,[name '.png']), 'Resolution',300);
close(fig);
end

function cs = controller_styles()
cs = {'C0_baseline',        'C0 structural ref.', '-.', 1.2, [0.122 0.467 0.706];
      'C1_hdiag_beta',      'C1 $h_{\rm diag}$-$\beta$', '--', 1.4, [0.900 0.590 0.350];
      'C2_full',            'C2 full (ours)',     '-',  1.7, [0.173 0.627 0.173];
      'C10_known_direction','C10 known-dir.',     ':',  1.4, [0.839 0.153 0.157];
      'C11_known_projected','C11 projected oracle','--',1.2, [0.45 0.33 0.62]};
end
