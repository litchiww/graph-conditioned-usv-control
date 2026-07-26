function cfg = usv_config()

cfg.default_params  = @default_params;
cfg.get_params      = @get_params;
cfg.case_list       = @case_list;
cfg.scenario_layout = @scenario_layout;
cfg.physical_map    = @physical_map;
cfg.mc_config       = @mc_config;
end

function p = default_params()

p.T            = 60.0;
p.dt           = 0.005;
p.save_stride  = 0.005;
p.solver       = 'ode15s';
p.rel_tol      = 1e-9;
p.abs_tol      = 1e-11;
p.max_step     = 0.005;
p.topology     = 'chain';
p.scenario     = 'offshore_nominal';
p.barrier_type = 'artanh';
p.integrate_barrier_state = false;
p.q_mode       = 'beta';
p.beta_mode    = 'opt';
p.beta_custom  = [2.1; 1.9; 1.0];
p.beta_rho_u   = 0.02;
p.beta_min     = 0.90;
p.beta_max     = 2.95;
p.beta_bar     = 5/3;
p.beta_step    = 0.05;
p.known_direction   = false;
p.command_limit     = Inf;
p.nonlinear_damping = true;
p.c1 = 2;       p.c2 = 32;

p.sigma1 = 0.05; p.sigma2 = 0.05;
p.search_deadzone1 = 0.0; p.search_deadzone2 = 0.0;
p.search_deadzone_mode = 'search_only';
p.search_deadzone_ramp_time = 0.0;
p.search_deadzone_power = 3;
p.search_lock_time = Inf;
p.search_lock_transition = 1.0;
p.search_lock_time1 = NaN; p.search_lock_time2 = NaN;
p.search_lock_transition1 = NaN; p.search_lock_transition2 = NaN;
p.search_lock_gain_min = 0.5;
p.search_lock_gain_full = 1.0;
p.search_lock_gain_min1 = NaN; p.search_lock_gain_full1 = NaN;
p.search_lock_gain_min2 = NaN; p.search_lock_gain_full2 = NaN;

p.search_lock_error_low1 = Inf; p.search_lock_error_high1 = Inf;
p.search_lock_error_low2 = Inf; p.search_lock_error_high2 = Inf;
p.search_lock_monitor_signal2 = 'chi';
p.search_lock_monitor_gain_low = 0.5;
p.search_lock_monitor_gain_high = 1.0;
p.search_lock_monitor_gain_low1 = NaN; p.search_lock_monitor_gain_high1 = NaN;
p.search_lock_monitor_gain_low2 = NaN; p.search_lock_monitor_gain_high2 = NaN;
p.search_sigma1_late = NaN;
p.search_sigma1_time = Inf;
p.search_sigma1_transition = 1.0;
p.search_c2_late = NaN;
p.search_c2_time = Inf;
p.search_c2_transition = 1.0;
p.search_command_bound1 = Inf; p.search_command_bound2 = Inf;
p.alpha_dot_regressor_mode = 'bound';
p.alpha_dot_regressor_blend = 0.25;

p.alpha_dot_filter_tau = 0.02;
p.alpha_dot_filter_margin = 0.0;
p.virtual_command_guard = false;
p.virtual_guard_ratio_low = 0.65;
p.virtual_guard_ratio_high = 0.80;
p.virtual_guard_contraction_width = 0.05;
p.second_layer_safety_flip = false;
p.state_derivative_filter = true;
p.state_derivative_filter_tau = 0.005;
p.direction_estimator = true;
p.direction_switch_mode = 'smooth';
p.direction_confidence1_min = 0.001;
p.direction_confidence1_full = 0.008;
p.direction_confidence2_min = 0.005;
p.direction_confidence2_full = 0.020;
p.safety_flip_ratio_low = 0.75;
p.safety_flip_ratio_high = 0.88;
p.safety_flip_trend_low = 0.20;
p.safety_flip_trend_high = 1.00;
p.safety_flip_contraction_low = 0.20;
p.safety_flip_contraction_high = 2.00;
p.safety_flip_release_low = 0.35;
p.safety_flip_release_high = 0.50;
p.safety_flip_on_rate = 50.0;
p.safety_flip_off_rate = 10.0;
p.envelope_filter_hold = false;
p.envelope_hold_eps = 1e-4;
p.envelope_filter_tau = 0.10;
p.envelope_decay_rate = Inf;
p.envelope_hold_activation = 0.05;
p.envelope_hold_end_time = Inf;
p.envelope_hold_gap_max = Inf;
p.gamma1 = 0.01; p.gamma2 = 0.01;
p.theta_leak1 = 0.005; p.theta_leak2 = 0.005;
p.rho = 1.0;    p.r = 0.49;
p.envelope_eps_abs = 1e-4;
p.envelope_eps_max = 1e-3;
p.eps1 = 0.20;  p.eps2 = 0.15;
p.nussbaum_family = 'xia2024';
p.nussbaum_alpha = 3.0;
p.nussbaum_beta = 0.50;
p.nussbaum_gamma = 0.0;
p.nussbaum_normalize_amplitude = false;
p.nussbaum_base = 2.0;
p.nussbaum_index = [1 2 3; 4 5 6];

p.qiao_a = ones(6,1);
p.qiao_b = 3*ones(6,1);
p.qiao_T = [0.04;0.06;0.08;0.010;0.015;0.020];

p.ma_amp = 1.0;
p.ma_alpha = 1.005;
p.ma_beta = 0.02;

p.yu_a = 0.9;
p.yu_b = 0.001;
p.yu_c = 0.7;
p.yu_omega = 0.4;
p.yu_phase = 'cos';
search_init=(-2.5:1:2.5)'*0.374;
p.search_init1 = search_init(1:3);
p.search_init2 = search_init(4:6);
p.theta_init1 = zeros(3,1);
p.theta_init2 = zeros(3,1);
p.lam0 = 0.05;  p.lamd = 0.15; p.lam_floor = 5e-4;
p.exp_clip = 45; p.ratio_clip = 1-1e-12;
p.barrier_margin = 1.0;

p.max_abs_state = 1e4; p.max_abs_u = 1e4;
p.max_ratio_guard = Inf;
p.g1_upper = [1.3; 1.3; 0.9];
p.g2_upper = [2.1; 2.1; 1.6];
p.proof_mu = 0.7;
p.young_eta = [4; 4; 4];
p.qa = 0.55; p.qb = 2.0; p.eps_q = 0.04;
p.kfa1 = 0.15; p.kfb1 = 0.04;
p.kfa2 = 0.02; p.kfb2 = 0.005;
p.probe_gain1 = 0.0; p.probe_gain2 = 0.0; p.probe_eps = 0.05;
p.mc_seed = 0;
p.mc_mode = 'none';
end

function names = case_list()

names = {'C0_baseline', 'C1_hdiag_beta', 'C2_full', ...
         'C3_star',     'C4_complete',     'C5_directed', ...
         'C6_rough_sea','C7_no_chi',       'C8_raw_error', ...
         'C9_decaying_diagnostic', 'C10_known_direction', 'C11_known_projected', ...
         'B1_no_damping'};
end

function p = get_params(name)

p = default_params();
switch name
    case 'C0_baseline'
        p.q_mode = 'hdiag'; p.beta_mode = 'hdiag';
        p.nonlinear_damping = false;
    case 'C1_hdiag_beta'
        p.q_mode = 'hdiag'; p.beta_mode = 'hdiag';
    case 'C2_full'

    case 'C3_star'
        p.topology = 'star';
    case 'C4_complete'
        p.topology = 'complete';
    case 'C5_directed'
        p.topology = 'directed';
    case 'C6_rough_sea'
        p.scenario = 'offshore_rough';
    case 'C7_no_chi'
        p.barrier_type = 'none';
    case 'C8_raw_error'
        p.q_mode = 'raw'; p.beta_mode = 'custom'; p.beta_custom = [1;1;1];
    case {'C9_decaying_diagnostic','C9_asymptotic'}
        p.scenario = 'decaying_diagnostic'; p.T = 150;
        p.lam_floor = 0;
        p.gamma1 = 0.001; p.gamma2 = 0.001;
        p.solver = 'rk4'; p.audit_detail = false; p.save_stride = 0.01;
    case 'C10_known_direction'
        p.known_direction = true;
        p.direction_estimator = false; p.state_derivative_filter = false;
    case 'C11_known_projected'
        p.known_direction = true;
        p.direction_estimator = false; p.state_derivative_filter = false;
        p.command_limit = 20.0;
    case 'B1_no_damping'
        p.nonlinear_damping = false;
    otherwise
        error('Unknown case: %s', name);
end
end

function s = scenario_layout()

s.lane_centers_y = [-800.0, 0.0, 800.0];
s.turbine_cols_x = [360.0, 1260.0, 2160.0];
s.turbine_offset_y_from_lane = 200.0;
s.safety_radius  = 80.0;
s.usv_min_separation = 200.0;
s.length_scale   = 400.0;
s.time_scale     = 12.0;
s.inspection_speed = 3.0;
s.path_progress_rate = s.time_scale * s.inspection_speed;
s.x_corridor     = [-100, 2400];
s.y_corridor     = [-1150, 1200];

turbines = [];
for col = s.turbine_cols_x
    for lane = s.lane_centers_y
        turbines = [turbines; col, lane + s.turbine_offset_y_from_lane];
    end
end
s.turbines = turbines;
end

function [X, Y, psi] = physical_map(t, x1_norm, x2_norm, row_id, layout)

Y_off = x1_norm * layout.length_scale;
Y = layout.lane_centers_y(row_id) + Y_off;
X = layout.path_progress_rate * t;
psi = atan2(x2_norm, 1.0);
end

function m = mc_config()

m.n_runs             = 30;
m.base_case          = 'C2_full';
m.groups = {'in_domain','sign_only','gain_only','initial_only','load_only','combined_stress'};
m.seed_offsets = [1000, 3000, 4000, 5000, 6000, 2000];
end
