function [out, arr, hp] = usv_simulate_3dof(case_dir)

if nargin < 1 || isempty(case_dir)
    here = fileparts(mfilename('fullpath'));
    case_dir = fullfile(fileparts(here), 'cases');
end

cfg = usv_config();
layout = cfg.scenario_layout();
ref = local_load_or_make_reference(case_dir);
[arr, hp] = local_run_plant(ref.arr, layout);
out = local_metrics(arr, hp);
end

function ref = local_load_or_make_reference(case_dir)
f = fullfile(case_dir, 'C2_full.mat');
if exist(f, 'file')
    ref = load(f);
    return;
end
cfg = usv_config();
p = cfg.get_params('C2_full');
[out, arr] = usv_simulate(p);
if ~exist(case_dir, 'dir'), mkdir(case_dir); end
save(f, 'out', 'arr', 'p', '-v7');
ref = struct('out', out, 'arr', arr, 'p', p);
end

function [arr, hp] = local_run_plant(ref, layout)
N = 3;
hp.N = N;
hp.time_scale = layout.time_scale;
hp.T = ref.t(end) * hp.time_scale;
hp.dt = 0.01;
hp.save_stride = 0.2;
hp.Ls = layout.length_scale;
hp.lane_y = layout.lane_centers_y(:);
hp.Uref = layout.inspection_speed * ones(N,1);
hp.psi_scale = 0.75;
hp.guidance_filter_s = 10.0;
hp.safety_sep = layout.usv_min_separation;
hp.control_rate_hz = 20.0;
hp.control_steps = round(1/(hp.control_rate_hz*hp.dt));
hp.packet_loss_rate = 0.05;
hp.random_seed = 314159;
hp.noise_y_m = 0.35;
hp.noise_psi_deg = 0.20;
hp.noise_r_degps = 0.15;
hp.noise_u_mps = 0.03;
hp.delta_max_deg = 25.0;

hp.m11 = [120; 126; 116];
hp.m22 = [185; 192; 178];
hp.m33 = [76; 82; 72];
hp.du  = [36; 38; 35];   hp.duu = [22; 24; 21];
hp.dv  = [92; 96; 88];   hp.dvv = [120; 126; 112];
hp.dr  = [28; 31; 26];   hp.drr = [48; 52; 45];
hp.Xmax = [820; 840; 800];
hp.Nmax = [235; 245; 225];
hp.Ku = [230; 240; 220];
hp.Ky = [0.030; 0.032; 0.028];
hp.Kpsi = [520; 545; 500];
hp.Kr = [135; 142; 128];
hp.wave_force = [10; 12; 9];
hp.wave_moment = [7; 8; 6];

refp = local_reference_pack(ref, layout, hp);

dt = hp.dt;
nsteps = round(hp.T/dt);
save_every = max(1, round(hp.save_stride/dt));
n_save = floor(nsteps/save_every) + 1;

X0 = [-80; 0; 80];
Y0 = refp.Y(1,:)';
psi0 = refp.psi(1,:)';
u0 = hp.Uref;
v0 = zeros(N,1);
r0 = refp.r(1,:)';
y = [X0; Y0; psi0; u0; v0; r0];

arr = local_alloc(n_save, N);
k_save = 0;
t = 0;
rng(hp.random_seed,'twister');
ctrl = local_init_controller(refp, hp);
for i = 0:nsteps
    if mod(i,hp.control_steps) == 0
        ctrl = local_update_controller(t, y, refp, hp, ctrl);
    end
    aux = local_aux(t, y, refp, hp, ctrl);
    if mod(i, save_every) == 0 || i == nsteps
        k_save = k_save + 1;
        arr = local_store(arr, k_save, t, y, aux);
    end
    if i == nsteps, break; end

    k1 = local_rhs(t,        y,             hp, ctrl.tauX, ctrl.tauN);
    k2 = local_rhs(t+dt/2.0, y + dt*k1/2.0, hp, ctrl.tauX, ctrl.tauN);
    k3 = local_rhs(t+dt/2.0, y + dt*k2/2.0, hp, ctrl.tauX, ctrl.tauN);
    k4 = local_rhs(t+dt,     y + dt*k3,     hp, ctrl.tauX, ctrl.tauN);
    y = y + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0;
    y(2*N+1:3*N) = atan2(sin(y(2*N+1:3*N)), cos(y(2*N+1:3*N)));
    t = t + dt;
end

arr = local_trim(arr, k_save);
arr.control_update_count = ctrl.update_count;
arr.yaw_saturation_count = ctrl.yaw_saturation_count;
arr.thrust_saturation_count = ctrl.thrust_saturation_count;
arr.max_abs_tauN_control = ctrl.max_abs_tauN;
arr.max_abs_tauX_control = ctrl.max_abs_tauX;
end

function refp = local_reference_pack(ref, layout, hp)
t0 = ref.t(:);
T = t0 <= hp.T / hp.time_scale + 1e-12;
t = t0(T) * hp.time_scale;
refp.t = t;
refp.X = hp.Uref(:)' .* t;
refp.Y = zeros(numel(t), hp.N);
win = max(1, round(hp.guidance_filter_s / mean(diff(t))));
refp.x1 = movmean(ref.x1(T,:), win, 1, 'Endpoints','shrink');

for j = 1:hp.N
    refp.Y(:,j) = layout.lane_centers_y(j) + hp.Ls * refp.x1(:,j);
end
refp.dY = local_gradient(refp.Y, t);
refp.psi = atan2(refp.dY, hp.Uref(:)');
refp.psi = max(min(refp.psi, 0.95*hp.psi_scale), -0.95*hp.psi_scale);
refp.r = local_gradient(refp.psi, t);
end

function ctrl = local_init_controller(refp, hp)
ctrl.Yref = refp.Y(1,:)';
ctrl.psiref = refp.psi(1,:)';
ctrl.rref = refp.r(1,:)';
ctrl.tauX = zeros(hp.N,1);
ctrl.tauN = zeros(hp.N,1);
ctrl.packet_ok = true(hp.N,1);
ctrl.packet_attempts = zeros(hp.N,1);
ctrl.packet_delivered = zeros(hp.N,1);
ctrl.update_count = 0;
ctrl.yaw_saturation_count = zeros(hp.N,1);
ctrl.thrust_saturation_count = zeros(hp.N,1);
ctrl.max_abs_tauN = zeros(hp.N,1);
ctrl.max_abs_tauX = zeros(hp.N,1);
end

function ctrl = local_update_controller(t, y, refp, hp, ctrl)
N = hp.N;
Y = y(N+1:2*N);
psi = y(2*N+1:3*N);
uu = y(3*N+1:4*N);
rr = y(5*N+1:6*N);

Ynew = local_interp(refp.t, refp.Y, t)';
psinew = local_interp(refp.t, refp.psi, t)';
rnew = local_interp(refp.t, refp.r, t)';
ctrl.update_count = ctrl.update_count + 1;
for j = 1:N
    ctrl.packet_attempts(j) = ctrl.packet_attempts(j) + 1;
    delivered = ctrl.update_count == 1 || rand >= hp.packet_loss_rate;
    ctrl.packet_ok(j) = delivered;
    if delivered
        ctrl.Yref(j) = Ynew(j);
        ctrl.psiref(j) = psinew(j);
        ctrl.rref(j) = rnew(j);
        ctrl.packet_delivered(j) = ctrl.packet_delivered(j) + 1;
    end

    Ym = Y(j) + hp.noise_y_m*randn;
    psim = psi(j) + deg2rad(hp.noise_psi_deg)*randn;
    rm = rr(j) + deg2rad(hp.noise_r_degps)*randn;
    um = uu(j) + hp.noise_u_mps*randn;

    eY = Ym - ctrl.Yref(j);
    psi_los = ctrl.psiref(j) - atan2(hp.Ky(j)*eY, hp.Uref(j));
    psi_los = local_sat(psi_los, 0.95*hp.psi_scale);
    epsi = atan2(sin(psi_los-psim), cos(psi_los-psim));
    er = ctrl.rref(j) - rm;
    Xff = hp.du(j)*hp.Uref(j) + hp.duu(j)*abs(hp.Uref(j))*hp.Uref(j);
    tauX_raw = Xff + hp.Ku(j)*(hp.Uref(j)-um);
    tauN_raw = hp.Kpsi(j)*epsi + hp.Kr(j)*er;
    ctrl.thrust_saturation_count(j) = ctrl.thrust_saturation_count(j) ...
        + double(abs(tauX_raw) >= hp.Xmax(j));
    ctrl.yaw_saturation_count(j) = ctrl.yaw_saturation_count(j) ...
        + double(abs(tauN_raw) >= hp.Nmax(j));
    ctrl.tauX(j) = local_sat(tauX_raw, hp.Xmax(j));
    ctrl.tauN(j) = local_sat(tauN_raw, hp.Nmax(j));
    ctrl.max_abs_tauX(j) = max(ctrl.max_abs_tauX(j), abs(ctrl.tauX(j)));
    ctrl.max_abs_tauN(j) = max(ctrl.max_abs_tauN(j), abs(ctrl.tauN(j)));
end
end

function dy = local_rhs(t, y, hp, tauX, tauN)
N = hp.N;
psi = y(2*N+1:3*N);
uu = y(3*N+1:4*N); vv = y(4*N+1:5*N); rr = y(5*N+1:6*N);

dX = zeros(N,1); dY = zeros(N,1); dpsi = zeros(N,1);
du = zeros(N,1); dv = zeros(N,1); dr = zeros(N,1);
for j = 1:N
    cj = cos(psi(j)); sj = sin(psi(j));
    dX(j) = cj*uu(j) - sj*vv(j);
    dY(j) = sj*uu(j) + cj*vv(j);
    dpsi(j) = rr(j);

    env = local_environment(t, j, hp);
    m11 = hp.m11(j); m22 = hp.m22(j); m33 = hp.m33(j);
    D = [hp.du(j) + hp.duu(j)*abs(uu(j));
         hp.dv(j) + hp.dvv(j)*abs(vv(j));
         hp.dr(j) + hp.drr(j)*abs(rr(j))];

    rhs1 = tauX(j) + env(1) + m22*vv(j)*rr(j) - D(1)*uu(j);
    rhs2 =           env(2) - m11*uu(j)*rr(j) - D(2)*vv(j);
    rhs3 = tauN(j) + env(3) + (m11-m22)*uu(j)*vv(j) - D(3)*rr(j);
    du(j) = rhs1 / m11;
    dv(j) = rhs2 / m22;
    dr(j) = rhs3 / m33;
end
dy = [dX; dY; dpsi; du; dv; dr];
end

function aux = local_aux(t, y, refp, hp, ctrl)
N = hp.N;
X = y(1:N); Y = y(N+1:2*N); psi = y(2*N+1:3*N);
uu = y(3*N+1:4*N); rr = y(5*N+1:6*N);

Yref = local_interp(refp.t, refp.Y, t);
psiref = local_interp(refp.t, refp.psi, t)';
rref = local_interp(refp.t, refp.r, t)';
Xref = hp.Uref(:) * t;

lat_err = Y - Yref(:);
Ncmd = ctrl.tauN;
Xcmd = ctrl.tauX;
delta_deg = hp.delta_max_deg * Ncmd ./ hp.Nmax;
power_W = abs(Xcmd.*uu) + abs(Ncmd.*rr);

aux.Xref = Xref;
aux.Yref = Yref(:);
aux.psiref = psiref;
aux.rref = rref;
aux.lat_err = lat_err;
aux.Ncmd = Ncmd;
aux.Xcmd = Xcmd;
aux.delta_deg = delta_deg;
aux.power_W = power_W;
aux.packet_ok = double(ctrl.packet_ok);
aux.packet_delivery_ratio = ctrl.packet_delivered ./ max(ctrl.packet_attempts,1);
aux.surge_err = hp.Uref(:) - uu;
aux.min_sep = local_min_sep(X, Y);
end

function env = local_environment(t, j, hp)
env = [0.35*hp.wave_force(j)*sin(0.23*t + 0.7*j);
       hp.wave_force(j)*sin(0.17*t + 0.4*j) + 0.45*hp.wave_force(j)*sin(0.41*t);
       hp.wave_moment(j)*sin(0.19*t + 0.9*j)];
end

function arr = local_alloc(n_save, N)
fields = {'X','Y','psi','u_surge','v_sway','r','Xref','Yref','psiref','rref', ...
          'lat_err','Ncmd','Xcmd','delta_deg','power_W', ...
          'packet_ok','packet_delivery_ratio','surge_err'};
arr.t = zeros(n_save,1);
for i = 1:numel(fields)
    arr.(fields{i}) = zeros(n_save,N);
end
arr.min_sep = zeros(n_save,1);
end

function arr = local_store(arr, k, t, y, aux)
N = size(arr.X,2);
arr.t(k) = t;
arr.X(k,:) = y(1:N)';
arr.Y(k,:) = y(N+1:2*N)';
arr.psi(k,:) = y(2*N+1:3*N)';
arr.u_surge(k,:) = y(3*N+1:4*N)';
arr.v_sway(k,:) = y(4*N+1:5*N)';
arr.r(k,:) = y(5*N+1:6*N)';
arr.Xref(k,:) = aux.Xref';
arr.Yref(k,:) = aux.Yref';
arr.psiref(k,:) = aux.psiref';
arr.rref(k,:) = aux.rref';
arr.lat_err(k,:) = aux.lat_err';
arr.Ncmd(k,:) = aux.Ncmd';
arr.Xcmd(k,:) = aux.Xcmd';
arr.delta_deg(k,:) = aux.delta_deg';
arr.power_W(k,:) = aux.power_W';
arr.packet_ok(k,:) = aux.packet_ok';
arr.packet_delivery_ratio(k,:) = aux.packet_delivery_ratio';
arr.surge_err(k,:) = aux.surge_err';
arr.min_sep(k) = aux.min_sep;
end

function arr = local_trim(arr, k_save)
fn = fieldnames(arr);
for i = 1:numel(fn)
    v = arr.(fn{i});
    arr.(fn{i}) = v(1:k_save,:);
end
end

function out = local_metrics(arr, hp)
tmask = arr.t >= arr.t(end)/2;
lat = arr.lat_err(tmask,:);
heading_err = atan2(sin(arr.psiref-arr.psi),cos(arr.psiref-arr.psi));
out.success = all(isfinite(lat(:))) && all(isfinite(arr.psi(:))) ...
    && min(arr.min_sep) > hp.safety_sep;
out.rms_lateral_m = sqrt(mean(lat(:).^2));
out.max_abs_lateral_m = max(abs(arr.lat_err(:)));
out.rms_heading_deg = rad2deg(sqrt(mean(heading_err(tmask,:).^2,'all')));
out.max_abs_heading_deg = rad2deg(max(abs(heading_err(:))));
out.min_sep_m = min(arr.min_sep);
out.max_yaw_moment_ratio = max(arr.max_abs_tauN_control ./ hp.Nmax);
out.max_thrust_ratio = max(arr.max_abs_tauX_control ./ hp.Xmax);
out.rms_rudder_deg = sqrt(mean(arr.delta_deg(tmask,:).^2,'all'));
out.max_abs_rudder_deg = hp.delta_max_deg*out.max_yaw_moment_ratio;
control_channel_samples = arr.control_update_count*hp.N;
out.yaw_saturation_count = sum(arr.yaw_saturation_count);
out.thrust_saturation_count = sum(arr.thrust_saturation_count);
out.control_channel_samples = control_channel_samples;
out.yaw_saturation_duty_pct = 100*out.yaw_saturation_count/control_channel_samples;
out.thrust_saturation_duty_pct = 100*out.thrust_saturation_count/control_channel_samples;
out.mean_total_power_W = mean(sum(arr.power_W,2));
out.energy_kWh = trapz(arr.t,sum(arr.power_W,2))/3.6e6;
out.packet_delivery_ratio = mean(arr.packet_delivery_ratio(end,:));
out.mean_speed_mps = mean(arr.u_surge(:));
out.max_heading_deg = max(abs(arr.psi(:))) * 180/pi;
out.T_final = arr.t(end);
out.rms_window = 'second half of physical horizon';
out.note = ['Sampled 3-DOF engineering validation with actuator limits, ' ...
            'sensor noise, lossy reference updates, and 20-Hz saturation accounting.'];
end

function y = local_interp(tgrid, values, t)
if isvector(values), values = values(:); end
y = interp1(tgrid, values, min(max(t, tgrid(1)), tgrid(end)), 'linear');
end

function g = local_gradient(x, t)
g = zeros(size(x));
for j = 1:size(x,2)
    g(:,j) = gradient(x(:,j), t);
end
end

function z = local_sat(x, lim)
z = max(min(x, lim), -lim);
end

function dmin = local_min_sep(X, Y)
N = numel(X);
dmin = Inf;
for a = 1:N
    for b = a+1:N
        dmin = min(dmin, hypot(X(a)-X(b), Y(a)-Y(b)));
    end
end
end
