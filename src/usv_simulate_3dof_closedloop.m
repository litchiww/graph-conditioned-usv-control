function [out, arr, hp] = usv_simulate_3dof_closedloop(overrides)

if nargin < 1, overrides = struct(); end
cfg = usv_config();
layout = cfg.scenario_layout();
[hp,p,api] = local_setup(cfg,layout,overrides);
[arr,run] = local_run(hp,api);
out = local_metrics(arr,hp,run);
end

function [hp,p,api] = local_setup(cfg,layout,overrides)
N=3;
controller_overrides=struct();
if isfield(overrides,'controller')
    controller_overrides=overrides.controller;
    overrides=rmfield(overrides,'controller');
end
hp.N=N;

hp.time_scale=6.0;
hp.T_normalized=60.0;
hp.T=hp.time_scale*hp.T_normalized;
hp.dt=0.01;
hp.save_stride=0.20;
hp.control_rate_hz=20.0;
hp.control_steps=round(1/(hp.control_rate_hz*hp.dt));
hp.lane_y=layout.lane_centers_y(:);
hp.safety_sep=layout.usv_min_separation;
hp.Uref=layout.inspection_speed*ones(N,1);
hp.psi_scale=1.25;
hp.r_scale=0.25;
hp.guidance_gain=[0.030;0.032;0.028];
hp.guidance_limit=0.45;
hp.guidance_wn=0.20;
hp.guidance_zeta=0.90;
hp.random_seed=271828;
hp.noise_psi_deg=0.05;
hp.noise_r_degps=0.03;
hp.initial_lateral_offset=[18;-14;11];
hp.initial_x1=[0.12;-0.08;0.10];
hp.initial_x2=[0.15;-0.12;0.10];
hp.probe_amplitude=[8;8;8];
hp.probe_frequency_hz=[0.40;0.47;0.54];
hp.probe_ramp_s=0.50;
hp.probe_hold_s=8.0;
hp.probe_fade_s=2.0;
hp.direction_confidence2_min=0.20;
hp.direction_confidence2_full=0.21;
hp.controller_c1=4.0;
hp.controller_gamma1=0.001;
hp.controller_gamma2=0.001;
hp.controller_theta_leak1=0.02;
hp.controller_theta_leak2=0.02;

hp.m11=[120;126;116];
hp.m22=[185;192;178];
hp.m33=[76;82;72];
hp.du=[36;38;35];       hp.duu=[22;24;21];
hp.dv=[92;96;88];       hp.dvv=[120;126;112];
hp.dr=[28;31;26];       hp.drr=[48;52;45];
hp.Xmax=[820;840;800];
hp.Nmax=[235;245;225];
hp.Ku=[230;240;220];
hp.wave_force=[10;12;9];
hp.wave_moment=[7;8;6];
hp.environment_ramp_s=30.0;

hp.yaw_polarity=[1;-1;1];
hp.g1=hp.time_scale*hp.r_scale/hp.psi_scale*ones(N,1);
hp.g2_magnitude=[2.0;2.0;1.5];
hp.g2=hp.yaw_polarity.*hp.g2_magnitude;
hp.tau_scale=hp.g2_magnitude.*hp.m33.*hp.r_scale/hp.time_scale;

names=fieldnames(overrides);
for k=1:numel(names)
    if ~isfield(hp,names{k})
        error('Unknown closed-loop 3-DOF override: %s',names{k});
    end
    hp.(names{k})=overrides.(names{k});
end
hp.T=hp.time_scale*hp.T_normalized;
hp.control_steps=round(1/(hp.control_rate_hz*hp.dt));
hp.g1=hp.time_scale*hp.r_scale/hp.psi_scale*ones(N,1);
hp.g2=hp.yaw_polarity.*hp.g2_magnitude;
hp.tau_scale=hp.g2_magnitude.*hp.m33.*hp.r_scale/hp.time_scale;

p=cfg.get_params('C2_full');
p.scenario='physical_yaw_closedloop';
p.g1_override=hp.g1;
p.g2_override=hp.g2;
p.g1_upper=1.08*abs(hp.g1);
p.g2_upper=1.05*abs(hp.g2);
p.direction_confidence2_min=hp.direction_confidence2_min;
p.direction_confidence2_full=hp.direction_confidence2_full;
p.c1=hp.controller_c1;
p.gamma1=hp.controller_gamma1;
p.gamma2=hp.controller_gamma2;
p.theta_leak1=hp.controller_theta_leak1;
p.theta_leak2=hp.controller_theta_leak2;
p.T=hp.T_normalized;
names=fieldnames(controller_overrides);
for k=1:numel(names)
    if ~isfield(p,names{k})
        error('Unknown controller override: %s',names{k});
    end
    p.(names{k})=controller_overrides.(names{k});
end
api=usv_simulate(p,'controller_api');
end

function [arr,run] = local_run(hp,api)
N=hp.N;
dt=hp.dt;
nsteps=round(hp.T/dt);
save_every=max(1,round(hp.save_stride/dt));
n_save=floor(nsteps/save_every)+1;

X0=[-80;0;80];
Y0=hp.lane_y+hp.initial_lateral_offset;
psi_cmd0=local_guidance_command(Y0,hp);
psi_d0=psi_cmd0;
r_d0=zeros(N,1);
x10=hp.initial_x1;
x20=hp.initial_x2;
psi0=psi_d0+hp.psi_scale*x10;
r0=r_d0+hp.r_scale*x20;
u0=hp.Uref;
v0=zeros(N,1);

y=[X0;Y0;psi0;u0;v0;r0;psi_d0;r_d0];
[z,~]=api.init(x10,x20);
ctrl=local_init_control(N,z);
arr=local_alloc(n_save,N);
k_save=0;
t=0;
rng(hp.random_seed,'twister');

for i=0:nsteps
    if mod(i,hp.control_steps)==0
        [ctrl,z]=local_update_control(t,y,z,hp,api,ctrl);
    end
    aux=local_aux(t,y,hp,ctrl);
    if mod(i,save_every)==0 || i==nsteps
        k_save=k_save+1;
        arr=local_store(arr,k_save,t,y,aux,ctrl);
    end
    if i==nsteps, break; end

    k1=local_rhs(t,       y,          hp,ctrl.tauX,ctrl.tauN);
    k2=local_rhs(t+dt/2, y+dt*k1/2, hp,ctrl.tauX,ctrl.tauN);
    k3=local_rhs(t+dt/2, y+dt*k2/2, hp,ctrl.tauX,ctrl.tauN);
    k4=local_rhs(t+dt,   y+dt*k3,   hp,ctrl.tauX,ctrl.tauN);
    y=y+dt*(k1+2*k2+2*k3+k4)/6;
    y(2*N+1:3*N)=local_wrap(y(2*N+1:3*N));
    y(6*N+1:7*N)=local_wrap(y(6*N+1:7*N));
    t=t+dt;
end

arr=local_trim(arr,k_save);
run.control_update_count=ctrl.update_count;
run.yaw_saturation_count=ctrl.yaw_saturation_count;
run.thrust_saturation_count=ctrl.thrust_saturation_count;
run.max_abs_tauN=ctrl.max_abs_tauN;
run.max_abs_tauX=ctrl.max_abs_tauX;
run.max_abs_u=ctrl.max_abs_u;
run.final_z=z;
end

function ctrl = local_init_control(N,z)
ctrl.z=z;
ctrl.tauX=zeros(N,1);
ctrl.tauN=zeros(N,1);
ctrl.u_cmd=zeros(N,1);
ctrl.u_probe=zeros(N,1);
ctrl.u_applied=zeros(N,1);
ctrl.sig=struct();
ctrl.update_count=0;
ctrl.yaw_saturation_count=zeros(N,1);
ctrl.thrust_saturation_count=zeros(N,1);
ctrl.max_abs_tauN=zeros(N,1);
ctrl.max_abs_tauX=zeros(N,1);
ctrl.max_abs_u=zeros(N,1);
end

function [ctrl,znew] = local_update_control(t,y,z,hp,api,ctrl)
N=hp.N;
psi=y(2*N+1:3*N);
uu=y(3*N+1:4*N);
rr=y(5*N+1:6*N);
psi_d=y(6*N+1:7*N);
r_d=y(7*N+1:8*N);

psi_m=psi+deg2rad(hp.noise_psi_deg)*randn(N,1);
r_m=rr+deg2rad(hp.noise_r_degps)*randn(N,1);
x1=local_wrap(psi_m-psi_d)/hp.psi_scale;
x2=(r_m-r_d)/hp.r_scale;
tau=t/hp.time_scale;
dtau=(1/hp.control_rate_hz)/hp.time_scale;

[sig1,k1,probe1]=local_controller_stage(api,tau,t,x1,x2,z,hp);
[~,k2]=local_controller_stage(api,tau+dtau/2,t+0.5/hp.control_rate_hz, ...
    x1,x2,z+dtau*k1/2,hp);
[~,k3]=local_controller_stage(api,tau+dtau/2,t+0.5/hp.control_rate_hz, ...
    x1,x2,z+dtau*k2/2,hp);
[~,k4]=local_controller_stage(api,tau+dtau,t+1/hp.control_rate_hz, ...
    x1,x2,z+dtau*k3,hp);
znew=z+dtau*(k1+2*k2+2*k3+k4)/6;

Xff=hp.du.*hp.Uref+hp.duu.*abs(hp.Uref).*hp.Uref;
tauX_raw=Xff+hp.Ku.*(hp.Uref-uu);
u_total=sig1.u+probe1;
tauN_raw=hp.yaw_polarity.*hp.tau_scale.*u_total;
tauX=local_sat(tauX_raw,hp.Xmax);
tauN=local_sat(tauN_raw,hp.Nmax);

ctrl.z=znew;
ctrl.sig=sig1;
ctrl.u_cmd=sig1.u;
ctrl.u_probe=probe1;
ctrl.u_applied=tauN./(hp.yaw_polarity.*hp.tau_scale);
ctrl.tauX=tauX;
ctrl.tauN=tauN;
ctrl.update_count=ctrl.update_count+1;
ctrl.yaw_saturation_count=ctrl.yaw_saturation_count+double(abs(tauN_raw)>=hp.Nmax);
ctrl.thrust_saturation_count=ctrl.thrust_saturation_count+double(abs(tauX_raw)>=hp.Xmax);
ctrl.max_abs_tauN=max(ctrl.max_abs_tauN,abs(tauN));
ctrl.max_abs_tauX=max(ctrl.max_abs_tauX,abs(tauX));
ctrl.max_abs_u=max(ctrl.max_abs_u,abs(u_total));
end

function [sig,dz,probe] = local_controller_stage(api,tau,t,x1,x2,z,hp)
[preview,~]=api.eval(tau,x1,x2,z);
probe=(1-preview.direction_blend2).*local_probe(t,hp);
[sig,dz]=api.eval_probed(tau,x1,x2,z,probe);
end

function dy = local_rhs(t,y,hp,tauX,tauN)
N=hp.N;
Y=y(N+1:2*N);
psi=y(2*N+1:3*N);
uu=y(3*N+1:4*N);
vv=y(4*N+1:5*N);
rr=y(5*N+1:6*N);
psi_d=y(6*N+1:7*N);
r_d=y(7*N+1:8*N);

psi_cmd=local_guidance_command(Y,hp);
dpsi_d=r_d;
dr_d=hp.guidance_wn^2.*local_wrap(psi_cmd-psi_d) ...
    -2*hp.guidance_zeta*hp.guidance_wn.*r_d;
dX=zeros(N,1); dY=zeros(N,1); dpsi=rr;
du=zeros(N,1); dv=zeros(N,1); dr=zeros(N,1);
for j=1:N
    cj=cos(psi(j)); sj=sin(psi(j));
    dX(j)=cj*uu(j)-sj*vv(j);
    dY(j)=sj*uu(j)+cj*vv(j);
    env=local_environment(t,j,hp);
    D=[hp.du(j)+hp.duu(j)*abs(uu(j));
       hp.dv(j)+hp.dvv(j)*abs(vv(j));
       hp.dr(j)+hp.drr(j)*abs(rr(j))];
    rhs1=tauX(j)+env(1)+hp.m22(j)*vv(j)*rr(j)-D(1)*uu(j);
    rhs2=env(2)-hp.m11(j)*uu(j)*rr(j)-D(2)*vv(j);
    rhs3=tauN(j)+env(3)+(hp.m11(j)-hp.m22(j))*uu(j)*vv(j)-D(3)*rr(j);
    du(j)=rhs1/hp.m11(j);
    dv(j)=rhs2/hp.m22(j);
    dr(j)=rhs3/hp.m33(j);
end
dy=[dX;dY;dpsi;du;dv;dr;dpsi_d;dr_d];
end

function psi_cmd = local_guidance_command(Y,hp)
cross_track=Y-hp.lane_y;
psi_cmd=-atan2(hp.guidance_gain.*cross_track,hp.Uref);
psi_cmd=local_sat(psi_cmd,hp.guidance_limit);
end

function probe = local_probe(t,hp)
if t<hp.probe_ramp_s
    window=sin(0.5*pi*t/hp.probe_ramp_s)^2;
elseif t<=hp.probe_hold_s
    window=1;
elseif t<hp.probe_hold_s+hp.probe_fade_s
    window=cos(0.5*pi*(t-hp.probe_hold_s)/hp.probe_fade_s)^2;
else
    window=0;
end
phase=[0;2*pi/3;4*pi/3];
probe=window*hp.probe_amplitude.*sin(2*pi*hp.probe_frequency_hz*t+phase);
end

function aux = local_aux(t,y,hp,ctrl)
N=hp.N;
X=y(1:N);
Y=y(N+1:2*N);
psi=y(2*N+1:3*N);
uu=y(3*N+1:4*N);
rr=y(5*N+1:6*N);
psi_d=y(6*N+1:7*N);
r_d=y(7*N+1:8*N);
x1=local_wrap(psi-psi_d)/hp.psi_scale;
x2=(rr-r_d)/hp.r_scale;
if isempty(fieldnames(ctrl.sig))
    error('Controller must be initialized before storing closed-loop data.');
end
aux.psi_cmd=local_guidance_command(Y,hp);
aux.x1=x1;
aux.x2=x2;
aux.lat_err=Y-hp.lane_y;
aux.heading_err=local_wrap(psi-psi_d);
aux.yaw_rate_err=rr-r_d;
aux.Xref=hp.Uref.*t;
aux.Yref=hp.lane_y;
aux.power_W=abs(ctrl.tauX.*uu)+abs(ctrl.tauN.*rr);
aux.min_sep=local_min_sep(X,Y);
end

function arr = local_alloc(n_save,N)
fields={'X','Y','psi','u_surge','v_sway','r','psi_d','r_d','psi_cmd', ...
    'x1','x2','lat_err','heading_err','yaw_rate_err','Xref','Yref', ...
    'e','q','alpha1','F','chi','ratio','u_cmd','u_probe','u_applied','tauN','tauX', ...
    's1','s2','eta1','eta2','direction_blend1','direction_blend2','power_W'};
arr.t=zeros(n_save,1);
for k=1:numel(fields)
    arr.(fields{k})=zeros(n_save,N);
end
arr.min_sep=zeros(n_save,1);
end

function arr = local_store(arr,k,t,y,aux,ctrl)
N=size(arr.X,2);
arr.t(k)=t;
arr.X(k,:)=y(1:N)';
arr.Y(k,:)=y(N+1:2*N)';
arr.psi(k,:)=y(2*N+1:3*N)';
arr.u_surge(k,:)=y(3*N+1:4*N)';
arr.v_sway(k,:)=y(4*N+1:5*N)';
arr.r(k,:)=y(5*N+1:6*N)';
arr.psi_d(k,:)=y(6*N+1:7*N)';
arr.r_d(k,:)=y(7*N+1:8*N)';
arr.psi_cmd(k,:)=aux.psi_cmd';
arr.x1(k,:)=aux.x1';
arr.x2(k,:)=aux.x2';
arr.lat_err(k,:)=aux.lat_err';
arr.heading_err(k,:)=aux.heading_err';
arr.yaw_rate_err(k,:)=aux.yaw_rate_err';
arr.Xref(k,:)=aux.Xref';
arr.Yref(k,:)=aux.Yref';
arr.e(k,:)=ctrl.sig.e';
arr.q(k,:)=ctrl.sig.q';
arr.alpha1(k,:)=ctrl.sig.alpha1';
arr.F(k,:)=ctrl.sig.F';
arr.chi(k,:)=ctrl.sig.chi';
arr.ratio(k,:)=ctrl.sig.ratio';
arr.u_cmd(k,:)=ctrl.u_cmd';
arr.u_probe(k,:)=ctrl.u_probe';
arr.u_applied(k,:)=ctrl.u_applied';
arr.tauN(k,:)=ctrl.tauN';
arr.tauX(k,:)=ctrl.tauX';
arr.s1(k,:)=ctrl.z(1:N)';
arr.s2(k,:)=ctrl.z(N+1:2*N)';
arr.eta1(k,:)=ctrl.sig.eta1';
arr.eta2(k,:)=ctrl.sig.eta2';
arr.direction_blend1(k,:)=ctrl.sig.direction_blend1';
arr.direction_blend2(k,:)=ctrl.sig.direction_blend2';
arr.power_W(k,:)=aux.power_W';
arr.min_sep(k)=aux.min_sep;
end

function arr = local_trim(arr,k_save)
names=fieldnames(arr);
for k=1:numel(names)
    arr.(names{k})=arr.(names{k})(1:k_save,:);
end
end

function out = local_metrics(arr,hp,run)
steady=arr.t>=arr.t(end)/2;
control_channel_samples=run.control_update_count*hp.N;
out.success=all(isfinite(arr.x1(:))) && all(isfinite(arr.x2(:))) ...
    && min(arr.min_sep)>hp.safety_sep && max(abs(arr.ratio(:)))<1;
out.rms_lateral_m=sqrt(mean(arr.lat_err(steady,:).^2,'all'));
out.max_abs_lateral_m=max(abs(arr.lat_err(:)));
out.rms_heading_deg=rad2deg(sqrt(mean(arr.heading_err(steady,:).^2,'all')));
out.max_abs_heading_deg=rad2deg(max(abs(arr.heading_err(:))));
out.rms_yaw_rate_degps=rad2deg(sqrt(mean(arr.yaw_rate_err(steady,:).^2,'all')));
out.rms_cooperative_error=sqrt(mean(sum(arr.e(steady,:).^2,2)));
out.max_envelope_ratio=max(abs(arr.ratio(:)));
out.min_envelope_margin_degps=rad2deg(min((arr.F-abs(arr.x2))*hp.r_scale,[],'all'));
out.min_sep_m=min(arr.min_sep);
out.max_yaw_moment_ratio=max(run.max_abs_tauN./hp.Nmax);
out.max_thrust_ratio=max(run.max_abs_tauX./hp.Xmax);
out.max_abs_u=max(run.max_abs_u);
out.yaw_saturation_count=sum(run.yaw_saturation_count);
out.thrust_saturation_count=sum(run.thrust_saturation_count);
out.control_channel_samples=control_channel_samples;
out.yaw_saturation_duty_pct=100*out.yaw_saturation_count/control_channel_samples;
out.thrust_saturation_duty_pct=100*out.thrust_saturation_count/control_channel_samples;
out.mean_total_power_W=mean(sum(arr.power_W,2));
out.energy_kWh=trapz(arr.t,sum(arr.power_W,2))/3.6e6;
out.final_direction_blend1=arr.direction_blend1(end,:);
out.final_direction_blend2=arr.direction_blend2(end,:);
out.final_identified_g2_sign=sign(arr.eta2(end,:));
out.final_control_polarity=-sign(arr.eta2(end,:));
out.actual_g2_sign=sign(hp.g2)';
out.direction_lock_time_s=local_lock_times(arr.t,arr.direction_blend2);
out.max_abs_probe=max(abs(arr.u_probe(:)));
out.first_layer_mapping_error=max(abs(hp.time_scale*hp.r_scale/hp.psi_scale-hp.g1));
mapped_g2=hp.time_scale*(hp.yaw_polarity.*hp.tau_scale)./(hp.m33*hp.r_scale);
out.input_mapping_error=max(abs(mapped_g2-hp.g2));
out.engineering_success=out.success && out.max_envelope_ratio<0.90 ...
    && out.yaw_saturation_count==0 ...
    && isequal(out.final_identified_g2_sign,out.actual_g2_sign);
out.T_final=arr.t(end);
out.note=['The theorem controller directly drives the 3-DOF yaw moment; ' ...
    'the prior offline-reference hydrodynamic check is retained separately.'];
end

function lock_time = local_lock_times(t,blend)
lock_time=NaN(1,size(blend,2));
for j=1:size(blend,2)
    k=find(blend(:,j)>=1-1e-6,1);
    if ~isempty(k), lock_time(j)=t(k); end
end
end

function env = local_environment(t,j,hp)
ramp=1-exp(-t/hp.environment_ramp_s);
env=ramp*[0.35*hp.wave_force(j)*sin(0.23*t+0.7*j);
          hp.wave_force(j)*sin(0.17*t+0.4*j)+0.45*hp.wave_force(j)*sin(0.41*t);
          hp.wave_moment(j)*sin(0.19*t+0.9*j)];
end

function z = local_sat(x,lim)
z=max(min(x,lim),-lim);
end

function a = local_wrap(a)
a=atan2(sin(a),cos(a));
end

function dmin = local_min_sep(X,Y)
N=numel(X);
dmin=Inf;
for i=1:N
    for j=i+1:N
        dmin=min(dmin,hypot(X(i)-X(j),Y(i)-Y(j)));
    end
end
end
