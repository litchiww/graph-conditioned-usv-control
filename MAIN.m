function MAIN(varargin)

HERE = fileparts(mfilename('fullpath'));
addpath(fullfile(HERE,'src'));
case_dir = fullfile(HERE,'cases');
fig_dir  = fullfile(HERE,'figures');
data_dir = fullfile(HERE,'data');
for d = {case_dir, fig_dir, data_dir}, if ~exist(d{1},'dir'), mkdir(d{1}); end, end

OPT.mc_mode      = 'relabel';
OPT.traj_rms     = 'steady';
OPT.keep_overview_consensus = true;

if isempty(varargin)
    do = struct('math',true,'solver',true,'signs',true,'sensitivity',true, ...
        'cases',true,'long',true,'mc',false,'lambda',true,'hydro',true, ...
        'closedloop',true,'plots',true,'animate',true);
else
    do = struct('math',false,'solver',false,'signs',false,'sensitivity',false, ...
        'cases',false,'long',false,'mc',false,'lambda',false,'hydro',false, ...
        'closedloop',false,'plots',false,'animate',false);
    stage = lower(varargin{1});
    if strcmp(stage,'realusv'), stage = 'hydro'; end
    if any(strcmp(stage,{'hydro_closedloop','direct3dof'})), stage = 'closedloop'; end
    if strcmp(stage,'audit')
        do.math=true; do.solver=true; do.signs=true; do.sensitivity=true; do.long=true;
        stage='';
    end
    if ~isempty(stage)
    do.(stage) = true;
    end
    if strcmp(stage,'plots')

    end
end

cfg = usv_config();

if do.math
    fprintf('\n=== Gate audit: mathematical identities and gain certificate ===\n');
    usv_validation_suite('math',HERE);
end
if do.solver
    fprintf('\n=== Gate audit: solver and step convergence ===\n');
    usv_validation_suite('solver',HERE);
end
if do.signs
    fprintf('\n=== Gate audit: exhaustive 64 fixed-sign patterns ===\n');
    usv_validation_suite('signs',HERE);
end
if do.sensitivity
    fprintf('\n=== Gate audit: leakage and denominator sensitivity ===\n');
    usv_validation_suite('sensitivity',HERE);
end

if do.cases
    fprintf('\n=== Stage 1: canonical and ablation cases ===\n');
    names = cfg.case_list();
    for i = 1:numel(names)
        nm = names{i};
        fprintf('  [%2d/%2d] %-22s ', i, numel(names), nm);
        tic;
        p = cfg.get_params(nm);
        [out, arr] = usv_simulate(p);
        fprintf('done (%.1fs)  RMS=%.4f  envelope_ratio=%.3f\n', toc, out.rms_e_steady, out.max_ratio);
        save(fullfile(case_dir, [nm '.mat']), 'out','arr','p','-v7');
    end

    extra_ss = {'SS_SS1','offshore_ss1'; 'SS_SS3','offshore_ss3'};
    for k = 1:size(extra_ss,1)
        fprintf('  [extra] %-22s ', extra_ss{k,1});
        tic; p = cfg.get_params('C2_full'); p.scenario = extra_ss{k,2};
        [out, arr] = usv_simulate(p);
        fprintf('done (%.1fs)  RMS=%.4f\n', toc, out.rms_e_steady);
        save(fullfile(case_dir, [extra_ss{k,1} '.mat']), 'out','arr','p','-v7');
    end
end

if do.long
    fprintf('\n=== Stage 1b: 600-s nominal boundedness audit ===\n');
    p = cfg.get_params('C2_full'); p.T = 600;
    [out, arr] = usv_simulate(p);
    save(fullfile(case_dir,'C2_long_600.mat'),'out','arr','p','-v7');
    fid = fopen(fullfile(data_dir,'long_horizon_summary.csv'),'w');
    if fid >= 0
        fprintf(fid,'horizon_s,rms_e_steady,max_ratio,max_abs_u,max_theta1,max_theta2,final_theta1,final_theta2,success\n');
        fprintf(fid,'%.0f,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%d\n', ...
            p.T,out.rms_e_steady,out.max_ratio,out.max_abs_u, ...
            max(arr.theta1,[],'all'),max(arr.theta2,[],'all'), ...
            max(arr.theta1(end,:)),max(arr.theta2(end,:)),out.success);
        fclose(fid);
    end
    fprintf('  RMS=%.4f ratio=%.3f max(theta)=[%.3f %.3f]\n', ...
        out.rms_e_steady,out.max_ratio,max(arr.theta1,[],'all'),max(arr.theta2,[],'all'));
end

if do.lambda
    fprintf('\n=== Stage 2: cross-topology lambda_beta scan ===\n');
    tools = usv_metrics();
    rows = tools.lambda_scan(40, 15.0, 42);
    save(fullfile(data_dir,'lambda_scan.mat'),'rows','-v7');

    local_lambda_scan_figure(rows, fig_dir);
end

if do.mc
    fprintf('\n=== Stage 3: Monte Carlo robustness ===\n');
    mcfg = cfg.mc_config();
    tools = usv_metrics();
    layout = cfg.scenario_layout();
    template = struct('group','','seed',NaN,'rms_e_steady',NaN,'path_rms_m',NaN,'max_ratio',NaN, ...
                      'paper_ok',false,'safe_ok',false,'fail_reason','none','fail_time',NaN, ...
                      'max_s',NaN,'max_F',NaN,'max_alpha1',NaN,'max_a2',NaN,'max_u',NaN);
    groups = mcfg.groups;
    offsets = mcfg.seed_offsets;
    mc = repmat(template, numel(groups)*mcfg.n_runs, 1);
    ix = 0;
    for gi = 1:numel(groups)
        fprintf('  Group: %s\n', strrep(groups{gi},'_','-'));
        for i = 1:mcfg.n_runs
            ix = ix + 1;
            seed = offsets(gi) + i;
            p = cfg.get_params(mcfg.base_case);
            p.T = 60; p.mc_seed = seed; p.mc_mode = groups{gi};
            mc(ix).group = groups{gi}; mc(ix).seed = seed;
            try
                [out,arr] = usv_simulate(p);
            catch ME
                warning('MC seed %d failed: %s',seed,ME.message);
                mc(ix).fail_reason = 'exception'; continue;
            end
            tm = tools.task_metrics(arr, layout);
            mc(ix).rms_e_steady = out.rms_e_steady;
            mc(ix).path_rms_m   = tm.path_rms_m;
            mc(ix).max_ratio    = out.max_ratio;
            mc(ix).paper_ok     = out.tracking_success && out.safety_success;
            mc(ix).safe_ok      = out.safety_success;
            mc(ix).fail_reason  = out.reason;
            mc(ix).fail_time    = arr.t(end);
            mc(ix).max_s        = max(abs([arr.s1(:); arr.s2(:)]));
            mc(ix).max_F        = max(arr.F(:));
            mc(ix).max_alpha1   = max(abs(arr.alpha1(:)));
            rr = min(abs(arr.ratio(:)), 1-1e-12);
            mc(ix).max_a2       = max(1./(1-rr.^2));
            mc(ix).max_u        = max(abs(arr.u(:)));
            if mod(i,5)==0 || ~out.success
                fprintf('    %2d/%d done   reason=%s\n', i, mcfg.n_runs, out.reason);
            end
        end
    end
    mc_schema_version = 30;
    save(fullfile(data_dir,'monte_carlo.mat'),'mc','mc_schema_version','-v7');
    local_write_mc_summary(mc, fullfile(data_dir,'monte_carlo_summary.csv'));
end

if do.hydro
    fprintf('\n=== Stage 4: 3-DOF hydrodynamic validation ===\n');
    tic;
    [out, arr, hp] = usv_simulate_3dof(case_dir);
    save(fullfile(data_dir,'hydro_validation.mat'),'out','arr','hp','-v7');
    local_write_hydro_summary(out, hp, fullfile(data_dir,'hydro_validation_summary.csv'));
    fprintf('  done (%.1fs)  lateral RMS=%.3f m  heading RMS=%.3f deg  min sep=%.1f m\n', ...
        toc, out.rms_lateral_m, out.rms_heading_deg, out.min_sep_m);
    usv_plots_3dof(struct('out',out,'arr',arr,'hp',hp), fig_dir);
end

if do.closedloop
    fprintf('\n=== Stage 4b: direct 3-DOF theorem-controller closure ===\n');
    tic;
    [out, arr, hp] = usv_simulate_3dof_closedloop();
    save(fullfile(data_dir,'hydro_closedloop_validation.mat'),'out','arr','hp','-v7.3');
    local_write_closedloop_summary(out,hp, ...
        fullfile(data_dir,'hydro_closedloop_validation_summary.csv'));
    fprintf(['  done (%.1fs)  lateral RMS=%.3f m  heading RMS=%.3f deg  ' ...
        'envelope=%.3f  saturation=%d/%d\n'], ...
        toc,out.rms_lateral_m,out.rms_heading_deg,out.max_envelope_ratio, ...
        out.yaw_saturation_count,out.control_channel_samples);
    usv_plots_3dof_closedloop(struct('out',out,'arr',arr,'hp',hp),fig_dir);
end

if do.plots
    fprintf('\n=== Stage 5: figures ===\n');
    R = local_load_cases(case_dir, [cfg.case_list(), {'SS_SS1','SS_SS3'}]);
    fprintf('Generating F1-F6 (control-theory figures):\n');
    usv_plots_control(R, fig_dir, OPT);

    R.SS_SS2 = R.C2_full;
    R.SS_SS4 = R.C6_rough_sea;

    mc = [];
    if exist(fullfile(data_dir,'monte_carlo.mat'),'file')
        S = load(fullfile(data_dir,'monte_carlo.mat'));
        if isfield(S,'mc_schema_version') && S.mc_schema_version==30
            mc = S.mc;
        else
            fprintf('  skipping legacy Monte Carlo data without schema-30 tag\n');
        end
    end
    fprintf('Generating A1-A5 (application figures):\n');
    usv_plots_app(R, mc, fig_dir, OPT);
    fprintf('Generating schematics:\n');
    usv_plots_schematic(R, case_dir, fig_dir);
    fprintf('Generating fixed-sign redesign audit:\n');
    usv_plots_redesign(HERE, fig_dir);
    if exist(fullfile(data_dir,'hydro_validation.mat'),'file')
        fprintf('Generating 3-DOF validation figure:\n');
        hydro = load(fullfile(data_dir,'hydro_validation.mat'));
        usv_plots_3dof(hydro, fig_dir);
    end
    if exist(fullfile(data_dir,'hydro_closedloop_validation.mat'),'file')
        fprintf('Generating direct 3-DOF closed-loop figure:\n');
        hydro_closedloop=load(fullfile(data_dir,'hydro_closedloop_validation.mat'));
        usv_plots_3dof_closedloop(hydro_closedloop,fig_dir);
    end

    local_write_summary_table(R, data_dir);
    local_write_task_metrics_csv(R, data_dir, OPT);
    local_write_ablation_csv(R, data_dir);
end

if do.animate
    fprintf('\n=== Stage 6: animation/video ===\n');
    R = local_load_cases(case_dir, {'C0_baseline','C2_full'});
    usv_animate(R, fig_dir, struct('fps',20,'gif',true,'mp4',true));
end

fprintf('\n=== All done.  See ./figures/ and ./data/. ===\n');
end

function R = local_load_cases(case_dir, names)
R = struct();
for i = 1:numel(names)
    nm = names{i};
    f = fullfile(case_dir, [nm '.mat']);
    if exist(f,'file'), S = load(f); R.(nm).out = S.out; R.(nm).arr = S.arr; R.(nm).p = S.p; end
end
end

function local_lambda_scan_figure(rows, fig_dir)
if isempty(rows), return; end
S = usv_style(); S.apply();
topo_colors = struct('chain',[0.122 0.467 0.706],'star',[1.0 0.498 0.055], ...
                     'complete',[0.173 0.627 0.173],'directed',[0.839 0.153 0.157]);
topo_markers = struct('chain','o','star','s','complete','^','directed','d');
fig = S.figure('double',0.52);
ax = axes('Parent',fig); hold(ax,'on');
all_lam=[]; all_rms=[];
for fn = fieldnames(topo_colors)'
    topo = fn{1};
    sub = rows(strcmp({rows.topology}, topo));
    if isempty(sub), continue; end
    lam = [sub.lambda_beta]; rms = [sub.rms_e];
    scatter(ax, lam, rms, 60, topo_colors.(topo), topo_markers.(topo), ...
        'filled','MarkerEdgeColor','k','LineWidth',0.4, ...
        'DisplayName',sprintf('%s (n=%d)',topo,numel(sub)),'MarkerFaceAlpha',0.6);
    all_lam = [all_lam lam]; all_rms = [all_rms rms];
end
mask = all_rms > 0 & all_lam > 0;
if sum(mask) > 10
    coef = polyfit(log10(all_lam(mask)), log10(all_rms(mask)), 1);
    lam_plot = logspace(log10(min(all_lam)), log10(max(all_lam)), 60);
    plot(ax, lam_plot, 10.^(coef(1)*log10(lam_plot)+coef(2)), 'k--','LineWidth',1.5, ...
        'DisplayName', sprintf('log-log fit (slope=%.2f)', coef(1)));
    corr = corrcoef(log10(all_lam(mask)), log10(all_rms(mask)));
    text(ax,0.04,0.93,{sprintf('Pearson $\\rho$(log,log) = %.3f',corr(1,2)), ...
        sprintf('$n=%d$ successful runs',numel(all_lam))}, ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7.5);
end
set(ax,'XScale','log','YScale','log');
xlabel(ax,'graph-conditioning index $\lambda_\beta=\lambda_{\min}(M_\beta^\top M_\beta)$');
ylabel(ax,'steady-state RMS $\|e\|$');
title(ax,'Cross-topology $\lambda_\beta$ scan','FontWeight','normal');
legend(ax,'show','Location','southeast'); grid(ax,'on'); ax.GridLineStyle=':'; ax.GridAlpha=0.3;
S.export(fig, fig_dir, 'fig_lambda_scan');
fprintf('  fig_lambda_scan\n');
end

function local_write_summary_table(R, data_dir)
tools = usv_metrics();
cfg   = usv_config();
names = cfg.case_list();
rows = struct('case',{},'track',{},'safe',{},'rms_e',{},'mean_e',{},'max_ratio',{}, ...
              'min_margin',{},'max_u',{},'lambda_beta',{},'topology',{},'scenario',{});
for i = 1:numel(names)
    nm = names{i};
    if ~isfield(R, nm), continue; end
    o = R.(nm).out;
    rows(end+1).case = nm;
    rows(end).track = o.tracking_success;
    rows(end).safe  = o.safety_success;
    rows(end).rms_e = o.rms_e_steady;
    rows(end).mean_e = o.mean_e_steady;
    rows(end).max_ratio  = o.max_ratio;
    rows(end).min_margin = o.min_safety_margin;
    rows(end).max_u      = o.max_abs_u;
    rows(end).lambda_beta = o.lambda_beta;
    rows(end).topology   = o.topology;
    rows(end).scenario   = o.scenario;
end
tools.write_table_md (rows, fullfile(data_dir,'summary_table.md'));
tools.write_table_csv(rows, fullfile(data_dir,'summary_table.csv'));
local_write_robustness_table(R, names, fullfile(data_dir,'robustness_table.tex'));
fprintf('  wrote summary_table.md, summary_table.csv, and robustness_table.tex\n');
end

function local_write_task_metrics_csv(R, data_dir, OPT)

tools  = usv_metrics();
cfg    = usv_config(); layout = cfg.scenario_layout();
keys   = {'C0_baseline','C1_hdiag_beta','C2_full','C10_known_direction','C11_known_projected'};
nice   = {'C0 hdiag/no damping','C1 hdiag-beta','C2 full (ours)', ...
          'C10 performance oracle','C11 projected oracle'};
fid = fopen(fullfile(data_dir,'task_metrics_table.csv'),'w');
if fid < 0
    warning('Could not write task metrics table.');
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'case,path_rms_m,normalized_effort_per_m,min_separation_m,ctrl_smooth\n');
for k=1:numel(keys)
    if ~isfield(R,keys{k}), continue; end
    tm = tools.task_metrics(R.(keys{k}).arr, layout);
    if strcmpi(OPT.traj_rms,'full')
        tm.path_rms_m = sqrt(mean(R.(keys{k}).arr.x1(:).^2)) * layout.length_scale;
    end
    fprintf(fid,'%s,%.6g,%.6g,%.6g,%.6g\n', nice{k}, ...
        tm.path_rms_m, tm.effort_per_m, tm.min_separation_m, tm.ctrl_smooth);
end
fprintf('  wrote task_metrics_table.csv\n');
end

function local_write_ablation_csv(R, data_dir)

keys = {'C2_full','B1_no_damping','C1_hdiag_beta', ...
        'C8_raw_error','C7_no_chi','C10_known_direction','C11_known_projected'};
nice = {'B0 full (ours)','B1 no damping','B2 hdiag-\beta', ...
        'B3 raw-e','B4 no-\chi','B5 known direction','B6 projected oracle'};
fid = fopen(fullfile(data_dir,'ablation_table.csv'),'w');
if fid < 0
    warning('Could not write ablation table.');
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'case,rms_e_steady,max_ratio,min_margin,max_abs_u,track,envelope\n');
for k=1:numel(keys)
    if ~isfield(R,keys{k}), continue; end
    o = R.(keys{k}).out;
    tr = 'FAIL'; if o.tracking_success, tr = 'OK'; end
    sf = 'FAIL'; if o.safety_success,  sf = 'OK'; end
    fprintf(fid,'%s,%.6g,%.6g,%.6g,%.6g,%s,%s\n', nice{k}, ...
        o.rms_e_steady, o.max_ratio, o.min_safety_margin, o.max_abs_u, tr, sf);
end
fprintf('  wrote ablation_table.csv\n');
end

function local_write_robustness_table(R, names, out_file)
fid = fopen(out_file,'w');
if fid < 0
    warning('Could not write robustness table: %s', out_file);
    return;
end
cleanup = onCleanup(@() fclose(fid));
for k = 1:numel(names)
    nm = names{k};
    if ~isfield(R, nm), continue; end

    if strcmp(nm,'C9_decaying_diagnostic'), continue; end
    o = R.(nm).out;
    fprintf(fid,'%s & %s & %.4g & %.3f & %.3f & %.2f & %s/%s \\\\\n', ...
        local_case_name(nm), o.topology, o.rms_e_steady, o.max_ratio, ...
        o.min_safety_margin, o.max_abs_u, local_ok(o.tracking_success), local_ok(o.safety_success));
end
fprintf(fid,'\\bottomrule\n');
end

function s = local_ok(tf)
if tf
    s = 'OK';
else
    s = 'FAIL';
end
end

function s = local_case_name(nm)
switch nm
    case 'C0_baseline'
        s = 'C0 hdiag/no damping';
    case 'C1_hdiag_beta'
        s = 'C1 hdiag-$\beta$';
    case 'C2_full'
        s = 'C2 full (ours)';
    case 'C3_star'
        s = 'C3 star';
    case 'C4_complete'
        s = 'C4 complete';
    case 'C5_directed'
        s = 'C5 directed';
    case 'C6_rough_sea'
        s = 'C6 rough sea';
    case 'C7_no_chi'
        s = 'C7 no-$\chi$';
    case 'C8_raw_error'
        s = 'C8 raw-$e$';
    case {'C9_decaying_diagnostic','C9_asymptotic'}
        s = 'C9 decaying diagnostic';
    case 'C10_known_direction'
        s = 'C10 known direction';
    case 'C11_known_projected'
        s = 'C11 projected oracle';
    case 'B1_no_damping'
        s = 'B1 no damping';
    otherwise
        s = strrep(nm,'_','\_');
end
end

function local_write_mc_summary(mc, out_file)
fid = fopen(out_file,'w');
if fid < 0, warning('Could not write %s',out_file); return; end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'group,total,completed,track_envelope,envelope_marginal,envelope_termination,numerical_termination,median_path_rms_m,max_completed_ratio,max_s,max_F,max_alpha1,max_a2,max_u\n');
groups = unique({mc.group},'stable');
for g = groups
    x = mc(strcmp({mc.group},g{1}));
    reason = {x.fail_reason};
    env = strcmp(reason,'constraint_violation');
    num = strcmp(reason,'numerical_limit') | strcmp(reason,'exception');
    done = find(~env & ~num);
    ok = sum([x(done).paper_ok]);
    safe = sum([x(done).safe_ok]);
    if isempty(done)
        med_rms = NaN; max_ratio = NaN;
    else
        med_rms = median([x(done).path_rms_m]);
        max_ratio = max([x(done).max_ratio]);
    end
    fprintf(fid,'%s,%d,%d,%d,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n',g{1},numel(x),numel(done), ...
        ok,safe-ok,sum(env),sum(num),med_rms,max_ratio,max([x.max_s],[],'omitnan'), ...
        max([x.max_F],[],'omitnan'),max([x.max_alpha1],[],'omitnan'), ...
        max([x.max_a2],[],'omitnan'),max([x.max_u],[],'omitnan'));
end
fprintf('  wrote monte_carlo_summary.csv\n');
end

function local_write_hydro_summary(out, hp, out_file)
fid = fopen(out_file,'w');
if fid < 0, warning('Could not write %s',out_file); return; end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,['control_rate_hz,target_packet_loss,actual_packet_loss,' ...
    'guidance_filter_s,rms_window_start_s,rms_lateral_m,max_lateral_m,' ...
    'rms_heading_deg,max_heading_error_deg,min_separation_m,rms_rudder_deg,' ...
    'max_rudder_deg,yaw_saturation_count,control_channel_samples,' ...
    'yaw_saturation_duty_pct,mean_fleet_power_kw,energy_kwh\n']);
fprintf(fid,['%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,' ...
    '%.6g,%.6g,%.6g,%.0f,%.0f,%.6g,%.6g,%.6g\n'], ...
    hp.control_rate_hz,hp.packet_loss_rate,1-out.packet_delivery_ratio, ...
    hp.guidance_filter_s,out.T_final/2,out.rms_lateral_m,out.max_abs_lateral_m, ...
    out.rms_heading_deg,out.max_abs_heading_deg,out.min_sep_m, ...
    out.rms_rudder_deg,out.max_abs_rudder_deg, ...
    out.yaw_saturation_count,out.control_channel_samples, ...
    out.yaw_saturation_duty_pct,out.mean_total_power_W/1000,out.energy_kWh);
fprintf('  wrote hydro_validation_summary.csv\n');
end

function local_write_closedloop_summary(out,hp,out_file)
fid=fopen(out_file,'w');
if fid<0, warning('Could not write %s',out_file); return; end
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,['time_scale_s,psi_scale_rad,r_scale_radps,control_rate_hz,' ...
    'heading_noise_deg,yaw_rate_noise_degps,environment_ramp_s,' ...
    'rms_window_start_s,rms_lateral_m,rms_heading_deg,rms_yaw_rate_degps,' ...
    'rms_cooperative_error,max_envelope_ratio,min_envelope_margin_degps,' ...
    'max_yaw_moment_ratio,yaw_saturation_count,control_channel_samples,' ...
    'yaw_saturation_duty_pct,max_abs_total_u,lock_time_1_s,lock_time_2_s,' ...
    'lock_time_3_s,identified_sign_1,identified_sign_2,identified_sign_3,' ...
    'first_layer_mapping_error,input_mapping_error,engineering_success\n']);
fprintf(fid,['%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,' ...
    '%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.0f,%.0f,%.6g,%.6g,' ...
    '%.6g,%.6g,%.6g,%.0f,%.0f,%.0f,%.6g,%.6g,%d\n'], ...
    hp.time_scale,hp.psi_scale,hp.r_scale,hp.control_rate_hz, ...
    hp.noise_psi_deg,hp.noise_r_degps,hp.environment_ramp_s,out.T_final/2, ...
    out.rms_lateral_m,out.rms_heading_deg,out.rms_yaw_rate_degps, ...
    out.rms_cooperative_error,out.max_envelope_ratio, ...
    out.min_envelope_margin_degps,out.max_yaw_moment_ratio, ...
    out.yaw_saturation_count,out.control_channel_samples, ...
    out.yaw_saturation_duty_pct,out.max_abs_u,out.direction_lock_time_s, ...
    out.final_identified_g2_sign,out.first_layer_mapping_error, ...
    out.input_mapping_error,out.engineering_success);
fprintf('  wrote hydro_closedloop_validation_summary.csv\n');
end
