function p=usv_redesign_candidate(p)

p.nussbaum_family='xia2024';
p.sigma1=0.05; p.sigma2=0.05; p.eps1=0.20; p.eps2=0.15;
p.search_deadzone1=0.0; p.search_deadzone2=0.0;
p.search_deadzone_mode='search_only';
p.search_sigma1_late=NaN; p.search_sigma1_time=Inf;

search_init=(-2.5:1:2.5)'*0.374;
p.search_init1=search_init(1:3); p.search_init2=search_init(4:6);
p.r=0.49;
p.integrate_barrier_state=false;
p.alpha_dot_regressor_mode='bound';
p.virtual_command_guard=false;
p.virtual_guard_ratio_low=0.45;
p.virtual_guard_ratio_high=0.60;
p.virtual_guard_contraction_width=0.005;
p.second_layer_safety_flip=false;
p.state_derivative_filter=true;
p.state_derivative_filter_tau=0.005;
p.direction_estimator=true;
p.direction_switch_mode='smooth';
p.direction_confidence1_min=0.001;
p.direction_confidence1_full=0.008;
p.direction_confidence2_min=0.005;
p.direction_confidence2_full=0.020;
p.safety_flip_ratio_low=0.75;
p.safety_flip_ratio_high=0.88;
p.safety_flip_trend_low=0.20;
p.safety_flip_trend_high=1.00;
p.safety_flip_contraction_low=0.20;
p.safety_flip_contraction_high=2.00;
p.safety_flip_release_low=0.35;
p.safety_flip_release_high=0.50;
p.safety_flip_on_rate=50.0;
p.safety_flip_off_rate=10.0;
p.envelope_filter_hold=false;
p.barrier_margin=1.0;
p.max_ratio_guard=0.98;
p.max_abs_u=100;

p.search_lock_time=Inf; p.search_lock_transition=1.0;
p.search_lock_gain_min=0.10; p.search_lock_gain_full=0.16;
p.search_lock_error_low1=0.10; p.search_lock_error_high1=0.30;
p.search_lock_monitor_gain_low=1e6; p.search_lock_monitor_gain_high=1e6+1;
p.search_lock_time1=Inf; p.search_lock_transition1=1.0;
p.search_lock_time2=Inf; p.search_lock_transition2=0.02;
p.search_lock_gain_min2=0.05; p.search_lock_gain_full2=0.10;
p.search_lock_monitor_signal2='chi';
p.search_lock_error_low2=Inf; p.search_lock_error_high2=Inf;

end
