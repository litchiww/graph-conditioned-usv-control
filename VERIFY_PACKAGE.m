function VERIFY_PACKAGE
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'src'));
patterns=readtable(fullfile(root,'data','sign_redesign_patterns.csv'),'TextType','string');
randomized=readtable(fullfile(root,'data','sign_redesign_mc.csv'),'TextType','string');
long_runs=readtable(fullfile(root,'data','sign_redesign_long.csv'),'TextType','string');
hydro=readtable(fullfile(root,'data','hydro_validation_summary.csv'));
hydro_closed=readtable(fullfile(root,'data','hydro_closedloop_validation_summary.csv'));
gain=readtable(fullfile(root,'data','gain_certificate.csv'));
math_audit=readtable(fullfile(root,'data','math_audit.csv'));
hydro_data=load(fullfile(root,'data','hydro_validation.mat'));
hydro_closed_data=load(fullfile(root,'data','hydro_closedloop_validation.mat'));
assert(height(patterns)==64);
assert(all(logical(patterns.completed)));
assert(all(logical(patterns.accepted)));
assert(max(patterns.max_ratio)<0.857);
assert(height(randomized)==30);
assert(all(logical(randomized.completed)));
assert(all(logical(randomized.sign_match)));
assert(all(logical(randomized.accepted)));
assert(max(randomized.max_ratio)<0.919);
assert(height(long_runs)==5);
assert(all(logical(long_runs.completed)));
assert(all(logical(long_runs.sign_match)));
assert(all(logical(long_runs.accepted)));
assert(all(logical(math_audit.pass)));
cfg=usv_config();
p=usv_redesign_candidate(cfg.get_params('C2_full'));
beta=gain.beta;
expected_cq=p.proof_mu*p.c1-(beta.*gain.g1_upper).^2./(2*gain.young_weight);
expected_cchi=p.proof_mu*p.c2-gain.young_weight/2;
expected_nu=min([2*min(expected_cq)*gain.lambda_beta(1), ...
    2*min(expected_cchi),p.theta_leak1,p.theta_leak2]);
assert(height(gain)==3);
assert(all(abs(gain.mu-p.proof_mu)<1e-12));
assert(max(abs(gain.c_q-expected_cq))<1e-10);
assert(max(abs(gain.c_chi-expected_cchi))<1e-10);
assert(all(gain.c_q>0) && all(gain.c_chi>0));
assert(max(abs(gain.nu-expected_nu))<1e-12);
mask=hydro_data.arr.t>=hydro_data.arr.t(end)/2;
lat_rms=sqrt(mean(hydro_data.arr.lat_err(mask,:).^2,'all'));
heading_error=atan2(sin(hydro_data.arr.psiref-hydro_data.arr.psi), ...
    cos(hydro_data.arr.psiref-hydro_data.arr.psi));
heading_rms=rad2deg(sqrt(mean(heading_error(mask,:).^2,'all')));
yaw_count=sum(hydro_data.arr.yaw_saturation_count);
control_samples=hydro_data.arr.control_update_count*hydro_data.hp.N;
yaw_duty=100*yaw_count/control_samples;
assert(abs(hydro.guidance_filter_s-hydro_data.hp.guidance_filter_s)<1e-12);
assert(abs(hydro.rms_window_start_s-hydro_data.arr.t(end)/2)<1e-9);
assert(abs(hydro.rms_lateral_m-lat_rms)<1e-5);
assert(abs(hydro.rms_heading_deg-heading_rms)<1e-5);
assert(hydro.yaw_saturation_count==yaw_count);
assert(hydro.control_channel_samples==control_samples);
assert(abs(hydro.yaw_saturation_duty_pct-yaw_duty)<1e-5);
assert(height(hydro_closed)==1);
assert(logical(hydro_closed.engineering_success));
assert(logical(hydro_closed_data.out.engineering_success));
assert(abs(hydro_closed.rms_lateral_m-2.37461)<1e-5);
assert(abs(hydro_closed.rms_heading_deg-5.92920)<1e-5);
assert(abs(hydro_closed.rms_yaw_rate_degps-1.05237)<1e-5);
assert(abs(hydro_closed.rms_cooperative_error-0.0326717)<1e-7);
assert(abs(hydro_closed.max_envelope_ratio-0.863809)<1e-6);
assert(hydro_closed.min_envelope_margin_degps>1.88);
assert(hydro_closed.max_yaw_moment_ratio<0.442);
assert(hydro_closed.yaw_saturation_count==0);
assert(hydro_closed.control_channel_samples==21603);
assert(hydro_closed.yaw_saturation_duty_pct==0);
assert(all([hydro_closed.identified_sign_1,hydro_closed.identified_sign_2, ...
    hydro_closed.identified_sign_3]==[1,-1,1]));
assert(hydro_closed.first_layer_mapping_error<1e-12);
assert(hydro_closed.input_mapping_error<1e-12);
assert(abs(hydro_closed_data.out.rms_lateral_m-hydro_closed.rms_lateral_m)<1e-5);
assert(abs(hydro_closed_data.out.max_envelope_ratio-hydro_closed.max_envelope_ratio)<1e-6);
p.T=0.10;
p.audit_detail=false;
p.save_stride=0.01;
[out,~]=usv_simulate(p);
assert(out.success);
assert(out.T_final>=0.099);
fprintf('Package verification passed.\n');
fprintf('Fixed signs: %d/64, worst ratio %.6f.\n',sum(patterns.accepted),max(patterns.max_ratio));
fprintf('Random trials: %d/30, worst ratio %.6f.\n',sum(randomized.accepted),max(randomized.max_ratio));
fprintf('Long audits: %d/5.\n',sum(long_runs.accepted));
fprintf('Gain certificate: mu %.3f, min c_q %.6f, min c_chi %.6f.\n', ...
    p.proof_mu,min(gain.c_q),min(gain.c_chi));
fprintf('3-DOF second-half RMS: %.6f m lateral, %.6f deg heading.\n', ...
    hydro.rms_lateral_m,hydro.rms_heading_deg);
fprintf('3-DOF saturation: %d/%d control updates (%.6f%%).\n', ...
    hydro.yaw_saturation_count,hydro.control_channel_samples, ...
    hydro.yaw_saturation_duty_pct);
fprintf(['Direct 3-DOF yaw loop: %.6f m lateral, %.6f deg heading, ' ...
    'envelope %.6f, saturation %d/%d.\n'], ...
    hydro_closed.rms_lateral_m,hydro_closed.rms_heading_deg, ...
    hydro_closed.max_envelope_ratio,hydro_closed.yaw_saturation_count, ...
    hydro_closed.control_channel_samples);
end
