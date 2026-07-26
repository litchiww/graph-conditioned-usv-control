function tools = usv_metrics()

tools.moving_rms      = @moving_rms;
tools.task_metrics    = @task_metrics;
tools.lane_deviation  = @lane_deviation;
tools.inter_distance  = @inter_distance;
tools.lambda_scan     = @lambda_scan;
tools.write_table_md  = @write_table_md;
tools.write_table_csv = @write_table_csv;
end

function rms = moving_rms(x, t, window_s)

dt = median(diff(t));
n = max(2, round(window_s / dt));
cs = [0; cumsum(x.^2)];
rms = zeros(size(x));
nh = floor(n/2);
for i = 1:numel(x)
    lo = max(1, i-nh); hi = min(numel(x), i+nh);
    rms(i) = sqrt((cs(hi+1)-cs(lo))/(hi-lo+1));
end
end

function tm = task_metrics(arr, layout)

t = arr.t; mask = t >= t(end)/2;

tm.path_rms_m   = sqrt(mean(arr.x1(mask,:).^2, 'all')) * layout.length_scale;
tm.path_max_m   = max(abs(arr.x1(:)))            * layout.length_scale;

effort_total    = trapz(t, sum(arr.u.^2, 2));
dist            = layout.path_progress_rate * t(end) * 3;
tm.effort_total = effort_total;
tm.effort_per_m = effort_total / max(dist, 1);

du = diff(arr.u, 1, 1) ./ diff(t);
tm.ctrl_smooth  = sqrt(mean(du.^2, 'all'));

d = inter_distance(arr, layout);
tm.min_separation_m = min([d.d12; d.d13; d.d23]);
end

function ld = lane_deviation(arr, layout)

ld.t = arr.t;
for j = 1:3
    ld.dev(:,j) = abs(arr.x1(:,j) * layout.length_scale);
end
ld.rms = sqrt(mean(ld.dev(:).^2));
ld.max = max(ld.dev(:));
end

function d = inter_distance(arr, layout)

[X1,Y1,~] = phys(arr.t, arr.x1(:,1), arr.x2(:,1), 1, layout);
[X2,Y2,~] = phys(arr.t, arr.x1(:,2), arr.x2(:,2), 2, layout);
[X3,Y3,~] = phys(arr.t, arr.x1(:,3), arr.x2(:,3), 3, layout);
d.d12 = sqrt((X1-X2).^2 + (Y1-Y2).^2);
d.d13 = sqrt((X1-X3).^2 + (Y1-Y3).^2);
d.d23 = sqrt((X2-X3).^2 + (Y2-Y3).^2);
d.t = arr.t;
end

function [X,Y,psi] = phys(t,x1,x2,row,layout)
Y = layout.lane_centers_y(row) + x1*layout.length_scale;
X = layout.path_progress_rate * t;
psi = atan2(x2, 1.0);
end

function rows = lambda_scan(n_per_topo, T_scan, seed)

if nargin<1, n_per_topo = 40; end
if nargin<2, T_scan = 15.0;   end
if nargin<3, seed = 42;       end
cfg = usv_config();
topos = {'chain','star','complete','directed'};
rng(seed);
rows = struct('topology',{},'beta',{},'lambda_beta',{},'cond_M',{},...
              'rms_e',{},'max_ratio',{});
ix = 1;
for ti = 1:numel(topos)
    topo = topos{ti};
    p_probe = cfg.get_params('C2_full');
    p_probe.topology = topo;
    [~,~,~,hdiag] = local_get_H(topo);
    tgt = mean(hdiag);
    lo = 0.75*tgt; hi = 1.35*tgt;
    kept = 0;
    for i = 1:n_per_topo
        b = lo + (hi-lo)*rand(3,1);
        b = b * (tgt/mean(b));
        b = max(min(b, 1.6*tgt), 0.5*tgt);
        p_try = p_probe;
        p_try.T = T_scan; p_try.dt = 0.01; p_try.save_stride = 0.1;
        p_try.beta_mode = 'custom';
        p_try.beta_custom = b;
        try
            [out,arr] = usv_simulate(p_try);
        catch
            continue;
        end
        if ~out.success || out.max_ratio >= 1.0, continue; end
        e_norm = vecnorm(arr.e,2,2);
        mask = arr.t >= T_scan/2;
        if sum(mask) < 2, continue; end
        rows(ix).topology    = topo;
        rows(ix).beta        = b;
        rows(ix).lambda_beta = out.lambda_beta;
        rows(ix).cond_M      = out.cond_M;
        rows(ix).rms_e       = sqrt(mean(e_norm(mask).^2));
        rows(ix).max_ratio   = out.max_ratio;
        ix = ix + 1; kept = kept + 1;
    end
    fprintf('  scan %-10s kept %d / %d\n', topo, kept, n_per_topo);
end
end

function [A,b,H,hdiag] = local_get_H(topo)
switch lower(topo)
    case 'chain',    A=[0 1 0;1 0 1;0 1 0]; b=[1;0;0];
    case 'star',     A=[0 1 1;1 0 0;1 0 0]; b=[1;0;0];
    case 'complete', A=ones(3)-eye(3);      b=[1;0;0];
    case 'directed', A=[0 0 1;1 0 0;0 1 0]; b=[1;0;0];
end
H = diag(sum(A,2)+b) - A; hdiag = diag(H);
end

function write_table_md(rows, filepath)

fid = fopen(filepath,'w');
fprintf(fid, '| Case | Track | Envelope | RMS_e (steady) | mean |e| | max |x2|/F | min margin | max |u| | lambda_beta |\n');
fprintf(fid, '|---|:-:|:-:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:numel(rows)
    r = rows(i);
    fprintf(fid, '| %s | %s | %s | %.4g | %.4g | %.4g | %.4g | %.4g | %.4g |\n', ...
        r.case, tf(r.track), tf(r.safe), r.rms_e, r.mean_e, r.max_ratio, ...
        r.min_margin, r.max_u, r.lambda_beta);
end
fclose(fid);
end

function s = tf(b)
if b, s='OK'; else, s='FAIL'; end
end

function write_table_csv(rows, filepath)
fid = fopen(filepath,'w');
fprintf(fid, 'case,track,envelope,rms_e_steady,mean_e_steady,max_ratio,min_margin,max_abs_u,lambda_beta,topology,scenario\n');
for i = 1:numel(rows)
    r = rows(i);
    fprintf(fid, '%s,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s,%s\n', ...
        r.case, r.track, r.safe, r.rms_e, r.mean_e, r.max_ratio, ...
        r.min_margin, r.max_u, r.lambda_beta, r.topology, r.scenario);
end
fclose(fid);
end
