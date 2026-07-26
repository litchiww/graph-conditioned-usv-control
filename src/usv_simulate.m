function [out, arr] = usv_simulate(p, mode)

N = 3;
if nargin < 2 || isempty(mode), mode = 'simulate'; end
if ~isfield(p,'solver'), p.solver = 'ode15s'; end
if ~isfield(p,'rel_tol'), p.rel_tol = 1e-9; end
if ~isfield(p,'abs_tol'), p.abs_tol = 1e-11; end
if ~isfield(p,'max_step'), p.max_step = p.dt; end
if ~isfield(p,'save_stride'), p.save_stride = p.max_step; end
if ~isfield(p,'audit_detail'), p.audit_detail = true; end
if ~isfield(p,'barrier_margin'), p.barrier_margin = 1.0; end
if ~(isscalar(p.barrier_margin) && isfinite(p.barrier_margin) ...
        && p.barrier_margin>0 && p.barrier_margin<=1)
    error('barrier_margin must be a finite scalar in (0,1].');
end
if ~isfield(p,'virtual_command_guard'), p.virtual_command_guard=false; end
if ~isfield(p,'state_derivative_filter'), p.state_derivative_filter=false; end
if ~isfield(p,'direction_estimator'), p.direction_estimator=false; end
if p.virtual_command_guard
    if p.integrate_barrier_state
        error('virtual_command_guard requires physical-state integration.');
    end
    if ~strcmpi(p.alpha_dot_regressor_mode,'command_filter')
        error('virtual_command_guard requires alpha_dot_regressor_mode=command_filter.');
    end
    if ~(0<=p.virtual_guard_ratio_low ...
            && p.virtual_guard_ratio_low<p.virtual_guard_ratio_high ...
            && p.virtual_guard_ratio_high<1)
        error('Virtual-command guard ratios must satisfy 0 <= low < high < 1.');
    end
    if ~(isscalar(p.virtual_guard_contraction_width) ...
            && p.virtual_guard_contraction_width>0)
        error('virtual_guard_contraction_width must be positive.');
    end
end
if p.state_derivative_filter
    if p.integrate_barrier_state
        error('state_derivative_filter requires physical-state integration.');
    end
    if ~(isscalar(p.state_derivative_filter_tau) ...
            && p.state_derivative_filter_tau>0)
        error('state_derivative_filter_tau must be positive.');
    end
end
if p.direction_estimator && ~p.state_derivative_filter
    error('direction_estimator requires state_derivative_filter.');
end
if p.direction_estimator
    if ~(0<=p.direction_confidence1_min ...
            && p.direction_confidence1_min<p.direction_confidence1_full ...
            && 0<=p.direction_confidence2_min ...
            && p.direction_confidence2_min<p.direction_confidence2_full)
        error('Direction-estimator confidence bands must satisfy 0 <= min < full.');
    end
    if ~any(strcmpi(p.direction_switch_mode,{'smooth','relay'}))
        error('direction_switch_mode must be smooth or relay.');
    end
end
if p.second_layer_safety_flip && ~p.state_derivative_filter
    error('second_layer_safety_flip requires state_derivative_filter.');
end
if p.second_layer_safety_flip
    if ~(0<=p.safety_flip_release_low ...
            && p.safety_flip_release_low<p.safety_flip_release_high ...
            && p.safety_flip_release_high<p.safety_flip_ratio_low)
        error('Safety-flip release ratios must lie below the trigger band.');
    end
    if ~(p.safety_flip_on_rate>0 && p.safety_flip_off_rate>0)
        error('Safety-flip transition rates must be positive.');
    end
end

[~, b, H, hdiag] = local_graph_mats(p.topology);
[beta, M] = local_conditioner(p, H, hdiag);
ref_map = H \ b;
sc = local_scenario(p.scenario, p);

if ~isfield(p,'search_init1'), p.search_init1=0.02*ones(N,1); end
if ~isfield(p,'search_init2'), p.search_init2=0.02*ones(N,1); end
if ~isfield(p,'theta_init1'), p.theta_init1=zeros(N,1); end
if ~isfield(p,'theta_init2'), p.theta_init2=zeros(N,1); end
if strcmpi(mode,'controller_api')
    if p.integrate_barrier_state
        error('controller_api requires physical-state integration.');
    end
    out = local_controller_api(p,N,M,beta,b,H,ref_map,sc);
    arr = struct();
    return;
elseif ~strcmpi(mode,'simulate')
    error('Unknown usv_simulate mode: %s',mode);
end
y0 = [sc.x10; sc.x20; p.search_init1(:); p.search_init2(:); ...
    p.theta_init1(:); p.theta_init2(:)];
if p.integrate_barrier_state && ~strcmp(p.barrier_type,'none')

    p0=p; p0.integrate_barrier_state=false; p0.envelope_filter_hold=false;
    p0.alpha_dot_regressor_mode='bound';
    sig0=local_control(0,y0,p0,N,M,beta,b,H,ref_map,sc);
    y0(N+1:2*N)=sig0.chi;
end
if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
    if ~(isscalar(p.alpha_dot_filter_tau) && p.alpha_dot_filter_tau>0)
        error('alpha_dot_filter_tau must be a positive scalar.');
    end

    p0=p; p0.alpha_dot_regressor_mode='bound'; p0.envelope_filter_hold=false;
    p0.state_derivative_filter=false;
    p0.second_layer_safety_flip=false;
    p0.direction_estimator=false;
    sig0=local_control(0,y0,p0,N,M,beta,b,H,ref_map,sc);
    y0=[y0;sig0.alpha1];
end
if p.state_derivative_filter
    y0=[y0;sc.x20];
end
if p.direction_estimator
    y0=[y0;sc.x10;zeros(2*N,1)];
end
if p.second_layer_safety_flip
    y0=[y0;zeros(N,1)];
end
if p.envelope_filter_hold
    if ~(isscalar(p.envelope_filter_tau) && p.envelope_filter_tau>0)
        error('envelope_filter_tau must be a positive scalar.');
    end
    if ~isfield(p,'envelope_decay_rate'), p.envelope_decay_rate=Inf; end
    if ~(isscalar(p.envelope_decay_rate) && p.envelope_decay_rate>0)
        error('envelope_decay_rate must be a positive scalar.');
    end
    p0=p; p0.envelope_filter_hold=false;
    sig0=local_control(0,y0,p0,N,M,beta,b,H,ref_map,sc);
    y0=[y0;sig0.F];
end
rhs = @(t,y) local_rhs(t,y,p,N,M,beta,b,H,ref_map,sc);
event_fun = @(t,y) local_events(t,y,p,N,M,beta,b,H,ref_map,sc);

switch lower(p.solver)
    case {'ode89','adaptive'}
        opts = odeset('RelTol',p.rel_tol,'AbsTol',p.abs_tol, ...
            'MaxStep',p.max_step,'Events',event_fun);
        if isfield(p,'max_wall_time') && isfinite(p.max_wall_time)
            opts=odeset(opts,'OutputFcn',@(tt,yy,flag) local_wall_output(tt,yy,flag,p.max_wall_time));
        end
        tspan = local_output_grid(p.T, p.save_stride);
        [t,Y,te,ye,ie] = ode89(rhs,tspan,y0,opts);
    case {'ode15s','radau'}
        opts = odeset('RelTol',p.rel_tol,'AbsTol',p.abs_tol, ...
            'MaxStep',p.max_step,'Events',event_fun);
        if isfield(p,'max_wall_time') && isfinite(p.max_wall_time)
            opts=odeset(opts,'OutputFcn',@(tt,yy,flag) local_wall_output(tt,yy,flag,p.max_wall_time));
        end
        tspan = local_output_grid(p.T, p.save_stride);
        [t,Y,te,ye,ie] = ode15s(rhs,tspan,y0,opts);
    case 'rk4'
        [t,Y,te,ye,ie] = local_rk4(rhs,event_fun,y0,p.T,p.dt,p.save_stride);
    otherwise
        error('Unknown solver: %s',p.solver);
end

if ~isempty(te) && (isempty(t) || abs(t(end)-te(end)) > 10*eps(max(1,te(end))))
    t = [t; te(end)];
    Y = [Y; ye(end,:)];
end

arr = local_collect(t,Y,p,N,M,beta,b,H,ref_map,sc);
out = local_metrics(arr,p,M,beta,sc,te,ie);
end

function api = local_controller_api(p,N,M,beta,b,H,ref_map,sc)
api.init = @(x1,x2) local_api_init(x1,x2,p,N,M,beta,b,H,ref_map,sc);
api.eval = @(t,x1,x2,z) local_api_eval(t,x1,x2,z,p,N,M,beta,b,H,ref_map,sc);
api.eval_probed = @(t,x1,x2,z,u_probe) ...
    local_api_eval(t,x1,x2,z,p,N,M,beta,b,H,ref_map,sc,u_probe);
api.N = N;
api.beta = beta;
api.M = M;
api.H = H;
api.b = b;
api.g1 = sc.g1;
api.g2 = sc.g2;
api.params = p;
end

function [z0,sig0] = local_api_init(x1,x2,p,N,M,beta,b,H,ref_map,sc)
x1=x1(:); x2=x2(:);
if numel(x1)~=N || numel(x2)~=N
    error('controller_api init expects N-element x1 and x2 vectors.');
end
y0=[x1;x2;p.search_init1(:);p.search_init2(:); ...
    p.theta_init1(:);p.theta_init2(:)];
if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
    p0=p; p0.alpha_dot_regressor_mode='bound'; p0.envelope_filter_hold=false;
    p0.state_derivative_filter=false;
    p0.second_layer_safety_flip=false;
    p0.direction_estimator=false;
    sig=local_control(0,y0,p0,N,M,beta,b,H,ref_map,sc);
    y0=[y0;sig.alpha1];
end
if p.state_derivative_filter
    y0=[y0;x2];
end
if p.direction_estimator
    y0=[y0;x1;zeros(2*N,1)];
end
if p.second_layer_safety_flip
    y0=[y0;zeros(N,1)];
end
if p.envelope_filter_hold
    p0=p; p0.envelope_filter_hold=false;
    sig=local_control(0,y0,p0,N,M,beta,b,H,ref_map,sc);
    y0=[y0;sig.F];
end
sig0=local_control(0,y0,p,N,M,beta,b,H,ref_map,sc);
z0=y0(2*N+1:end);
end

function [sig,dz] = local_api_eval(t,x1,x2,z,p,N,M,beta,b,H,ref_map,sc,u_probe)
y=[x1(:);x2(:);z(:)];
sig=local_control(t,y,p,N,M,beta,b,H,ref_map,sc);
if nargin>=13 && ~isempty(u_probe) && p.direction_estimator
    sig.deta2=(1-sig.direction_blend2).*(sig.u+u_probe(:)).*sig.dx2_filter;
end
dz=[sig.ds1;sig.ds2;sig.dtheta1;sig.dtheta2];
if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
    dz=[dz;sig.dalpha_filter];
end
if p.state_derivative_filter
    dz=[dz;sig.dx2_filter];
end
if p.direction_estimator
    dz=[dz;sig.dx1_filter;sig.deta1;sig.deta2];
end
if p.second_layer_safety_flip
    dz=[dz;sig.dsafety_mode2];
end
if p.envelope_filter_hold
    dz=[dz;sig.dF_filter];
end
if numel(dz)~=numel(z)
    error('controller_api internal-state layout mismatch (%d versus %d).',numel(dz),numel(z));
end
end

function stop = local_wall_output(~,~,flag,limit_s)
persistent clock_token
stop=false;
if strcmp(flag,'init'), clock_token=tic; return; end
if strcmp(flag,'done'), clock_token=[]; return; end
if ~isempty(clock_token) && toc(clock_token)>limit_s, stop=true; end
end

function grid = local_output_grid(T, stride)
grid = (0:stride:T)';
if isempty(grid) || grid(end) < T-10*eps(max(1,T)), grid(end+1,1) = T; end
end

function [beta,M] = local_conditioner(p,H,hdiag)
N = size(H,1);
switch p.q_mode
    case 'raw'
        beta = ones(N,1); M = eye(N);
    case 'hdiag'
        beta = hdiag(:); M = diag(1./beta)*H';
    case 'beta'
        switch p.beta_mode
            case 'hdiag', beta = hdiag(:);
            case 'custom', beta = p.beta_custom(:);
            case 'opt', beta = local_opt_beta(H,hdiag,p);
            otherwise, error('Bad beta_mode: %s',p.beta_mode);
        end
        M = diag(1./beta)*H';
    otherwise
        error('Bad q_mode: %s',p.q_mode);
end
end

function dy = local_rhs(t,y,p,N,M,beta,b,H,ref_map,sc)
sig = local_control(t,y,p,N,M,beta,b,H,ref_map,sc);
x1 = y(1:N); x2 = sig.x2;
dx1 = sc.g1.*x2 + sc.f1(x1,x2,t);
dx2 = sc.g2.*sig.u + sc.f2(x1,x2,t);
if p.integrate_barrier_state && ~strcmp(p.barrier_type,'none')
    dsecond=sig.ax.*dx2-sig.aalpha.*sig.alpha1_dot+sig.Gamma.*sig.F_dot;
else
    dsecond=dx2;
end
dy = [dx1; dsecond; sig.ds1; sig.ds2; sig.dtheta1; sig.dtheta2];
if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
    dy=[dy;sig.dalpha_filter];
end
if p.state_derivative_filter
    dy=[dy;sig.dx2_filter];
end
if p.direction_estimator
    dy=[dy;sig.dx1_filter;sig.deta1;sig.deta2];
end
if p.second_layer_safety_flip
    dy=[dy;sig.dsafety_mode2];
end
if p.envelope_filter_hold
    dy=[dy;sig.dF_filter];
end
end

function [value,isterminal,direction] = local_events(t,y,p,N,M,beta,b,H,ref_map,sc)
sig = local_control(t,y,p,N,M,beta,b,H,ref_map,sc);
s = y(2*N+1:4*N);
if strcmp(p.barrier_type,'none') || p.integrate_barrier_state
    envelope_distance = 1;
else
    envelope_distance = min(sig.F-abs(sig.x2));
end
cap_distance = p.exp_clip-max(local_nussbaum_exponent(s,p));
state_distance = p.max_abs_state-max(abs(y));
if isfinite(p.max_abs_u), command_distance = p.max_abs_u-max(abs(sig.u));
else, command_distance = 1; end
if isfinite(p.max_ratio_guard)
    ratio_distance=p.max_ratio_guard-max(abs(sig.ratio));
else
    ratio_distance=1;
end
value = [envelope_distance; cap_distance; state_distance; command_distance; ratio_distance];
isterminal = ones(5,1);
direction = -ones(5,1);
end

function sig = local_control(t,y,p,N,M,beta,b,H,ref_map,sc)
x1 = y(1:N); second = y(N+1:2*N);
s1 = y(2*N+1:3*N); s2 = y(3*N+1:4*N);
th1 = y(4*N+1:5*N); th2 = y(5*N+1:6*N);
state_offset=6*N;
if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
    alpha_filter=y(state_offset+1:state_offset+N);
    state_offset=state_offset+N;
else
    alpha_filter=zeros(N,1);
end
if p.state_derivative_filter
    x2_filter=y(state_offset+1:state_offset+N);
    state_offset=state_offset+N;
else
    x2_filter=zeros(N,1);
end
if p.direction_estimator
    x1_filter=y(state_offset+1:state_offset+N);
    eta1=y(state_offset+N+1:state_offset+2*N);
    eta2=y(state_offset+2*N+1:state_offset+3*N);
    state_offset=state_offset+3*N;
else
    x1_filter=zeros(N,1); eta1=zeros(N,1); eta2=zeros(N,1);
end
if p.second_layer_safety_flip
    safety_mode2=max(0,min(1,y(state_offset+1:state_offset+N)));
    state_offset=state_offset+N;
else
    safety_mode2=zeros(N,1);
end
if p.envelope_filter_hold
    envelope_filter=y(state_offset+1:state_offset+N);
else
    envelope_filter=zeros(N,1);
end

yd = sc.yd(t); yd1 = sc.yd1(t); yd2 = sc.yd2(t);
e = H*x1-b*yd;
if strcmp(p.q_mode,'raw'), q=e; else, q=M*e; end

ell = p.lam_floor+p.lam0*exp(-p.lamd*t);
ell_dot_abs = p.lam0*p.lamd*exp(-p.lamd*t);
rho1 = 1+x1.^2;
if strcmp(p.q_mode,'raw')
    R1 = 1+abs(H)*rho1+abs(b)*abs(yd1);
else
    R1 = 1+beta.*rho1+abs(beta.*ref_map)*abs(yd1);
end

[damp1,damp1_prime] = local_damp(q,p,1);
base1 = p.c1*q+damp1;
z1 = q.*R1;
den1 = sqrt(z1.^2+ell^2+1e-18);
adapt1 = th1.*q.*R1.^2./den1;
bar1 = base1+adapt1;

[N1,Np1] = local_nussbaum(s1,p,1);
sigma1_eff=local_search_sigma(t,p,1,N);
P1raw=q.*bar1;
[deadzone_scale,deadzone_scale_dot]=local_deadzone_schedule(t,p);
d1=p.search_deadzone1*deadzone_scale;
d1_dot=p.search_deadzone1*deadzone_scale_dot;
d2=p.search_deadzone2*deadzone_scale;
matched_rational=strcmpi(p.search_deadzone_mode,'matched_rational');
matched_power=strcmpi(p.search_deadzone_mode,'matched_power');
matched_hard=strcmpi(p.search_deadzone_mode,'matched_hard');
matched_compact=strcmpi(p.search_deadzone_mode,'matched_compact');
matched_deadzone=matched_rational || matched_power || matched_hard || matched_compact;
if matched_rational && d1>0
    P1=P1raw;
    h1=P1./(P1+d1);
    dead_drive1=bar1.*h1;
elseif matched_power && d1>0
    power=p.search_deadzone_power;
    if ~(isscalar(power) && power>=3 && power==floor(power))
        error('search_deadzone_power must be an integer not smaller than 3.');
    end
    P1=P1raw; rpower=power-1;
    h1=(P1./(P1+d1)).^rpower;
    dead_drive1=bar1.*h1;
elseif matched_hard && d1>0
    P1=P1raw;
    m1=max(P1-d1,0);
    h1=zeros(N,1); active1=P1>d1;
    h1(active1)=m1(active1)./P1(active1);
    dead_drive1=bar1.*h1;
elseif matched_compact && d1>0
    P1=P1raw;
    m1=zeros(N,1); mp1=zeros(N,1);
    mid1=P1>d1 & P1<2*d1; high1=P1>=2*d1;
    m1(mid1)=(P1(mid1)-d1).^2./(2*d1);
    mp1(mid1)=(P1(mid1)-d1)./d1;
    m1(high1)=P1(high1)-1.5*d1; mp1(high1)=1;
    h1=zeros(N,1); active1=P1>d1;
    h1(active1)=m1(active1)./P1(active1);
    h1prime=zeros(N,1);
    h1prime(active1)=(mp1(active1).*P1(active1)-m1(active1))./P1(active1).^2;
    dead_drive1=bar1.*h1;
else
    P1=max(P1raw,0); h1=ones(N,1); dead_drive1=bar1;
end
[N1_effective,direction_blend1,direction_blend1_eta]= ...
    local_direction_override(N1,eta1,p,1);
[lock_blend1,lock_blend1_time_dot,lock_blend1_N,lock_blend1_q]= ...
    local_search_lock_blend(t,p,N1_effective,q,1);
locked_sign1=sign(N1_effective); locked_sign1(locked_sign1==0)=1;
identified_sign1=-sign(eta1); identified_sign1(identified_sign1==0)=1;
raw_search_drive1=N1.*dead_drive1;
raw_search1=(1-direction_blend1).*raw_search_drive1 ...
    +direction_blend1.*identified_sign1.*bar1;
raw1=(1-lock_blend1).*raw_search1+lock_blend1.*locked_sign1.*bar1;
scale1=sqrt(1+(raw1./p.search_command_bound1).^2);
search_scale1=sqrt(1+(raw_search_drive1./p.search_command_bound1).^2);
drive1=dead_drive1./search_scale1;
if p.known_direction
    alpha1 = -sign(sc.g1).*bar1./beta;
    ds1 = zeros(N,1);
else
    alpha1 = raw1./(scale1.*beta);
    if matched_deadzone
        ds1 = (1-direction_blend1).*(1-lock_blend1).*sigma1_eff.*q.*drive1;
    else
        ds1 = (1-direction_blend1).*(1-lock_blend1).*sigma1_eff.*max(q.*drive1-p.search_deadzone1,0);
    end
end
virtual_guard_blend=zeros(N,1);
if p.virtual_command_guard
    [Ftrial,~,~,~]=local_envelope(t,alpha1,zeros(N,1),p);
    [Fmemory,~,~,~]=local_envelope(t,alpha_filter,zeros(N,1),p);
    ratio_pressure=local_unit_smooth(abs(second)./Ftrial, ...
        p.virtual_guard_ratio_low,p.virtual_guard_ratio_high);
    contraction_pressure=local_unit_smooth(Fmemory-Ftrial,0, ...
        p.virtual_guard_contraction_width);
    virtual_guard_blend=ratio_pressure.*contraction_pressure;
    alpha1=(1-virtual_guard_blend).*alpha1+virtual_guard_blend.*alpha_filter;
end
dalpha_filter=zeros(N,1);
if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
    dalpha_filter=(alpha1-alpha_filter)./p.alpha_dot_filter_tau;
end
dtheta1 = p.gamma1*(z1.^2)./den1-p.theta_leak1*th1;
[F,Falpha,Ffilter,Fcmd] = local_envelope(t,alpha1,envelope_filter,p);
if p.envelope_filter_hold
    dF_filter=(Fcmd-envelope_filter)./p.envelope_filter_tau;
    dF_filter=max(dF_filter,-p.envelope_decay_rate);
else
    dF_filter=zeros(N,1);
end
pa = alpha1./F;
if p.integrate_barrier_state && ~strcmp(p.barrier_type,'none')
    [Ba,aalpha] = local_barrier(pa,p);
    Bx = second./F+Ba;
    ratio=local_barrier_inverse(Bx,p);
    x2=F.*ratio;
else
    x2=second;
    ratio=x2./F;
end
if p.state_derivative_filter
    dx2_filter=(x2-x2_filter)./p.state_derivative_filter_tau;
else
    dx2_filter=zeros(N,1);
end
if p.direction_estimator
    dx1_filter=(x1-x1_filter)./p.state_derivative_filter_tau;
else
    dx1_filter=zeros(N,1);
end
Fdot_est=Falpha.*dalpha_filter+Ffilter.*dF_filter;
ratio_trend_est=sign(x2).*dx2_filter./F-abs(x2).*Fdot_est./F.^2;

dx1=sc.g1.*x2+sc.f1(x1,x2,t);
if strcmp(p.q_mode,'raw')
    q_dot=H*dx1-b*yd1;
else
    q_dot=M*(H*dx1-b*yd1);
end
lock_blend1_dot=lock_blend1_time_dot+lock_blend1_N.*Np1.*ds1 ...
    +lock_blend1_q.*q_dot;
rho1_dot=2*x1.*dx1;
yd1_abs_dot=sign(yd1)*yd2;
if strcmp(p.q_mode,'raw')
    R1_dot=abs(H)*rho1_dot+abs(b)*yd1_abs_dot;
else
    R1_dot=beta.*rho1_dot+abs(beta.*ref_map)*yd1_abs_dot;
end
ell_dot=-p.lam0*p.lamd*exp(-p.lamd*t);
z1_dot=q_dot.*R1+q.*R1_dot;
hz1_exact=ell^2./den1.^3;
hell1_exact=-z1*ell./den1.^3;
h1=z1./den1;
adapt1_dot=dtheta1.*R1.*h1+th1.*R1_dot.*h1 ...
    +th1.*R1.*(hz1_exact.*z1_dot+hell1_exact*ell_dot);
bar1_dot=(p.c1+damp1_prime).*q_dot+adapt1_dot;
if p.known_direction
    alpha1_dot=-sign(sc.g1).*bar1_dot./beta;
else
    if matched_rational && d1>0
        P1_dot=q_dot.*bar1+q.*bar1_dot;
        h1_dot=(d1.*P1_dot-P1.*d1_dot)./(P1+d1).^2;
        dead_drive1_dot=bar1_dot.*h1+bar1.*h1_dot;
    elseif matched_power && d1>0
        P1_dot=q_dot.*bar1+q.*bar1_dot;
        h1_P=rpower*d1.*P1.^(rpower-1)./(P1+d1).^(rpower+1);
        h1_d=-rpower*P1.^rpower./(P1+d1).^(rpower+1);
        h1_dot=h1_P.*P1_dot+h1_d.*d1_dot;
        dead_drive1_dot=bar1_dot.*h1+bar1.*h1_dot;
    elseif matched_hard && d1>0
        P1_dot=q_dot.*bar1+q.*bar1_dot;
        h1_dot=zeros(N,1);
        h1_dot(active1)=(d1.*P1_dot(active1)-P1(active1).*d1_dot)./P1(active1).^2;
        dead_drive1_dot=bar1_dot.*h1+bar1.*h1_dot;
    elseif matched_compact && d1>0
        P1_dot=q_dot.*bar1+q.*bar1_dot;
        h1_d=zeros(N,1);
        h1_d(mid1)=(-P1(mid1).^2./(2*d1^2)+0.5)./P1(mid1);
        h1_d(high1)=-1.5./P1(high1);
        h1_dot=h1prime.*P1_dot+h1_d.*d1_dot;
        dead_drive1_dot=bar1_dot.*h1+bar1.*h1_dot;
    else
        dead_drive1_dot=bar1_dot;
    end
    raw_search_drive1_dot=Np1.*ds1.*dead_drive1+N1.*dead_drive1_dot;
    deta1_for_derivative=(1-direction_blend1).*x2.*dx1_filter;
    direction_blend1_dot=direction_blend1_eta.*deta1_for_derivative;
    raw_search1_dot=(1-direction_blend1).*raw_search_drive1_dot ...
        +direction_blend1.*identified_sign1.*bar1_dot ...
        +direction_blend1_dot.*(identified_sign1.*bar1-raw_search_drive1);
    raw1_dot=(1-lock_blend1).*raw_search1_dot ...
        +lock_blend1.*locked_sign1.*bar1_dot ...
        +lock_blend1_dot.*(locked_sign1.*bar1-raw_search1);
    alpha1_dot=raw1_dot./(beta.*(1+(raw1./p.search_command_bound1).^2).^(3/2));
end
F_dot=Falpha.*alpha1_dot+Ffilter.*dalpha_filter;

qdot_reg = 1+abs(M*H)*(p.g1_upper.*abs(x2)+rho1)+abs(M*b)*abs(yd1);
rho1dot_reg = 1+2*abs(x1).*(p.g1_upper.*abs(x2)+rho1);
if strcmp(p.q_mode,'raw')
    R1dot_reg = 1+abs(H)*rho1dot_reg+abs(b)*abs(yd2);
else
    R1dot_reg = 1+beta.*rho1dot_reg+abs(beta.*ref_map)*abs(yd2);
end

hz1=ell^2./den1.^3;
hell1=abs(z1*ell)./den1.^3;
zdot_reg=R1.*qdot_reg+abs(q).*R1dot_reg;
adapt_dot_reg=abs(dtheta1).*R1+abs(th1).*R1dot_reg ...
    +abs(th1).*R1.*(hz1.*zdot_reg+hell1*ell_dot_abs);
bar1dot_reg = (p.c1+abs(damp1_prime)).*qdot_reg+adapt_dot_reg;
if p.known_direction
    Areg = 1+bar1dot_reg./beta;
else
    if matched_rational && d1>0
        P1dot_reg=qdot_reg.*abs(bar1)+abs(q).*bar1dot_reg;
        h1dot_reg=(d1.*P1dot_reg+P1.*abs(d1_dot))./(P1+d1).^2;
        dead_drive1dot_reg=h1.*bar1dot_reg+abs(bar1).*h1dot_reg;
    elseif matched_power && d1>0
        P1dot_reg=qdot_reg.*abs(bar1)+abs(q).*bar1dot_reg;
        h1dot_reg=abs(h1_P).*P1dot_reg+abs(h1_d).*abs(d1_dot);
        dead_drive1dot_reg=h1.*bar1dot_reg+abs(bar1).*h1dot_reg;
    elseif matched_hard && d1>0
        P1dot_reg=qdot_reg.*abs(bar1)+abs(q).*bar1dot_reg;
        h1dot_reg=zeros(N,1);
        h1dot_reg(active1)=(d1.*P1dot_reg(active1)+P1(active1).*abs(d1_dot))./P1(active1).^2;
        dead_drive1dot_reg=h1.*bar1dot_reg+abs(bar1).*h1dot_reg;
    elseif matched_compact && d1>0
        P1dot_reg=qdot_reg.*abs(bar1)+abs(q).*bar1dot_reg;
        h1dot_reg=abs(h1prime).*P1dot_reg+abs(h1_d).*abs(d1_dot);
        dead_drive1dot_reg=h1.*bar1dot_reg+abs(bar1).*h1dot_reg;
    else
        dead_drive1dot_reg=bar1dot_reg;
    end
    raw_search_drive1dot_reg=abs(Np1.*ds1.*dead_drive1) ...
        +abs(N1).*dead_drive1dot_reg;
    deta1_for_reg=(1-direction_blend1).*x2.*dx1_filter;
    direction_blend1dot_reg=abs(direction_blend1_eta.*deta1_for_reg);
    Areg=1+((1-direction_blend1).*raw_search_drive1dot_reg ...
        +direction_blend1.*bar1dot_reg ...
        +direction_blend1dot_reg.*(abs(bar1)+abs(raw_search_drive1)))./beta;
end
if strcmpi(p.alpha_dot_regressor_mode,'exact_diagnostic')

    Areg=1+abs(alpha1_dot);
elseif strcmpi(p.alpha_dot_regressor_mode,'blended_diagnostic')
    exact_reg=1+abs(alpha1_dot);
    blend=p.alpha_dot_regressor_blend;
    Areg=exact_reg+blend.*max(Areg-exact_reg,0);
elseif strcmpi(p.alpha_dot_regressor_mode,'command_filter')

    Areg=1+abs(dalpha_filter)+p.alpha_dot_filter_margin;
else
    dalpha_filter=zeros(N,1);
end
if ~exist('dalpha_filter','var'), dalpha_filter=zeros(N,1); end

if strcmp(p.barrier_type,'none')
    chi = x2-alpha1;
    ax = ones(N,1); aalpha = ones(N,1); Gamma = zeros(N,1);
elseif p.integrate_barrier_state
    chi=second;
    [~,ax]=local_barrier(ratio,p);
    Gamma=Bx-Ba-ax.*ratio+aalpha.*pa;
else
    [Bx,ax] = local_barrier(x2./F,p);
    [Ba,aalpha] = local_barrier(pa,p);
    chi = F.*(Bx-Ba);
    Gamma = Bx-Ba-ax.*(x2./F)+aalpha.*pa;
end

Freg = abs(Falpha).*Areg+abs(Ffilter).*(1+abs(dalpha_filter));
rho2 = 1+x1.^2+x2.^2;
R2 = 1+ax.*rho2+aalpha.*Areg+abs(Gamma).*Freg;
[damp2,~] = local_damp(chi,p,2);
c2_eff=local_c2_schedule(t,p,N);
base2 = c2_eff.*chi+damp2;
z2 = chi.*R2;
den2 = sqrt(z2.^2+ell^2+1e-18);
adapt2 = th2.*chi.*R2.^2./den2;
bar2 = base2+adapt2;

[N2,Np2] = local_nussbaum(s2,p,2);
if p.second_layer_safety_flip
    ratio_pressure2=local_unit_smooth(abs(ratio), ...
        p.safety_flip_ratio_low,p.safety_flip_ratio_high);
    trend_pressure2=local_unit_smooth(ratio_trend_est, ...
        p.safety_flip_trend_low,p.safety_flip_trend_high);
    contraction_pressure2=local_unit_smooth(-Fdot_est, ...
        p.safety_flip_contraction_low,p.safety_flip_contraction_high);
    safety_trigger2=ratio_pressure2.*trend_pressure2.*(1-contraction_pressure2);
    safety_release2=1-local_unit_smooth(abs(ratio), ...
        p.safety_flip_release_low,p.safety_flip_release_high);
    dsafety_mode2=p.safety_flip_on_rate.*safety_trigger2.*(1-safety_mode2) ...
        -p.safety_flip_off_rate.*safety_release2.*safety_mode2;
    safety_flip2=safety_mode2;
else
    safety_flip2=zeros(N,1);
    safety_trigger2=zeros(N,1);
    dsafety_mode2=zeros(N,1);
end
N2_search=(1-2*safety_flip2).*N2;
[N2_effective,direction_blend2,~]=local_direction_override(N2_search,eta2,p,2);
sigma2_eff=p.sigma2.*ones(N,1);
P2raw=chi.*bar2;
if matched_rational && d2>0
    P2=P2raw;
    h2=P2./(P2+d2);
    dead_drive2=bar2.*h2;
elseif matched_power && d2>0
    power=p.search_deadzone_power;
    if ~(isscalar(power) && power>=3 && power==floor(power))
        error('search_deadzone_power must be an integer not smaller than 3.');
    end
    P2=P2raw; rpower=power-1;
    h2=(P2./(P2+d2)).^rpower;
    dead_drive2=bar2.*h2;
elseif matched_hard && d2>0
    P2=P2raw;
    m2=max(P2-d2,0);
    h2=zeros(N,1); active2=P2>d2;
    h2(active2)=m2(active2)./P2(active2);
    dead_drive2=bar2.*h2;
elseif matched_compact && d2>0
    P2=P2raw;
    m2=zeros(N,1);
    mid2=P2>d2 & P2<2*d2; high2=P2>=2*d2;
    m2(mid2)=(P2(mid2)-d2).^2./(2*d2);
    m2(high2)=P2(high2)-1.5*d2;
    h2=zeros(N,1); active2=P2>d2;
    h2(active2)=m2(active2)./P2(active2);
    dead_drive2=bar2.*h2;
else
    P2=max(P2raw,0); dead_drive2=bar2;
end
locked_sign2=sign(N2_effective); locked_sign2(locked_sign2==0)=1;
if strcmpi(p.search_lock_monitor_signal2,'ratio')
    lock_monitor2=ratio;
elseif strcmpi(p.search_lock_monitor_signal2,'chi')
    lock_monitor2=chi;
elseif strcmpi(p.search_lock_monitor_signal2,'safety')
    lock_monitor2=safety_mode2;
else
    error('search_lock_monitor_signal2 must be chi, ratio, or safety.');
end
[lock_blend2,~,~,~]=local_search_lock_blend(t,p,N2_effective,lock_monitor2,2);
identified_sign2=-sign(eta2); identified_sign2(identified_sign2==0)=1;
raw_search_drive2=N2_search.*dead_drive2;
raw_search2=(1-direction_blend2).*raw_search_drive2 ...
    +direction_blend2.*identified_sign2.*bar2;
raw2=(1-lock_blend2).*raw_search2+lock_blend2.*locked_sign2.*bar2;
scale2=sqrt(1+(raw2./p.search_command_bound2).^2);
search_scale2=sqrt(1+(raw_search_drive2./p.search_command_bound2).^2);
drive2=dead_drive2./search_scale2;
if p.known_direction
    u = -sign(sc.g2).*bar2./ax;
    ds2 = zeros(N,1);
else
    u = raw2./(scale2.*ax);
    if matched_deadzone
        ds2 = (1-direction_blend2).*(1-lock_blend2).*sigma2_eff.*chi.*drive2;
    else
        ds2 = (1-direction_blend2).*(1-lock_blend2).*sigma2_eff.*max(chi.*drive2-p.search_deadzone2,0);
    end
end
u = max(min(u,p.command_limit),-p.command_limit);
if p.direction_estimator
    deta1=(1-direction_blend1).*x2.*dx1_filter;
    deta2=(1-direction_blend2).*u.*dx2_filter;
else
    deta1=zeros(N,1); deta2=zeros(N,1);
end
dtheta2 = p.gamma2*(z2.^2)./den2-p.theta_leak2*th2;
phi2=ax.*sc.f2(x1,x2,t)-aalpha.*alpha1_dot+Gamma.*F_dot;

sig = struct('x2',x2,'e',e,'q',q,'alpha1',alpha1,'alpha1_dot',alpha1_dot, ...
    'F',F,'Falpha',Falpha,'Ffilter',Ffilter,'F_dot',F_dot,'phi2',phi2, ...
    'chi',chi,'u',u,'ratio',x2./F,'N1',N1,'N2',N2,'Np1',Np1,'Np2',Np2, ...
    'safety_flip2',safety_flip2,'N2_effective',N2_effective, ...
    'safety_trigger2',safety_trigger2,'dsafety_mode2',dsafety_mode2, ...
    'N1_effective',N1_effective,'direction_blend1',direction_blend1, ...
    'direction_blend2',direction_blend2, ...
    'R1',R1,'R2',R2,'ax',ax,'aalpha',aalpha,'Gamma',Gamma, ...
    'Areg',Areg,'Freg',Freg,'bar1',bar1,'bar2',bar2,'drive1',drive1,'drive2',drive2, ...
    'ds1',ds1,'ds2',ds2, ...
    'dtheta1',dtheta1,'dtheta2',dtheta2,'ell',ell, ...
    'alpha_filter',alpha_filter,'dalpha_filter',dalpha_filter, ...
    'x2_filter',x2_filter,'dx2_filter',dx2_filter, ...
    'x1_filter',x1_filter,'dx1_filter',dx1_filter, ...
    'eta1',eta1,'eta2',eta2,'deta1',deta1,'deta2',deta2, ...
    'ratio_trend_est',ratio_trend_est, ...
    'virtual_guard_blend',virtual_guard_blend, ...
    'envelope_filter',envelope_filter,'dF_filter',dF_filter, ...
    'lock_blend1',lock_blend1,'lock_blend2',lock_blend2, ...
    'sigma1_eff',sigma1_eff,'sigma2_eff',sigma2_eff);
sig.c2_eff=c2_eff;
end

function arr = local_collect(t,Y,p,N,M,beta,b,H,ref_map,sc)
n = numel(t);
arr.t=t; arr.x1=Y(:,1:N);
if p.integrate_barrier_state && ~strcmp(p.barrier_type,'none')
    arr.chi_state=Y(:,N+1:2*N); arr.x2=zeros(n,N);
else
    arr.x2=Y(:,N+1:2*N);
end
arr.s1=Y(:,2*N+1:3*N); arr.s2=Y(:,3*N+1:4*N);
arr.theta1=Y(:,4*N+1:5*N); arr.theta2=Y(:,5*N+1:6*N);
fieldsN = {'e','q','alpha1','F','chi','u','ratio','N1','N2','Np1','Np2', ...
    'R1','R2','ax','aalpha','Gamma','Areg','Freg','bar1','bar2','drive1','drive2','ds1','ds2', ...
    'dtheta1','dtheta2','f1','f2','edot','residual1', ...
    'alpha1_dot','F_dot','chi_dot_dir','alpha_filter','dalpha_filter', ...
    'virtual_guard_blend','safety_flip2','N2_effective', ...
    'safety_trigger2','dsafety_mode2','N1_effective', ...
    'direction_blend1','direction_blend2', ...
    'x2_filter','dx2_filter','x1_filter','dx1_filter','ratio_trend_est', ...
    'eta1','eta2','deta1','deta2', ...
    'envelope_filter','dF_filter', ...
    'lock_blend1','lock_blend2','sigma1_eff','sigma2_eff','c2_eff'};
for k=1:numel(fieldsN), arr.(fieldsN{k})=zeros(n,N); end
arr.yd=zeros(n,1); arr.ell=zeros(n,1); arr.graph_identity_error=zeros(n,1);

for k=1:n
    y=Y(k,:)'; tk=t(k);
    sig=local_control(tk,y,p,N,M,beta,b,H,ref_map,sc);
    x1=y(1:N); x2=sig.x2;
    arr.x2(k,:)=x2';
    f1=sc.f1(x1,x2,tk); f2=sc.f2(x1,x2,tk);
    edot=H*(sc.g1.*x2+f1)-b*sc.yd1(tk);
    if p.audit_detail
        dx2=sc.g2.*sig.u+f2;
        if p.integrate_barrier_state && ~strcmp(p.barrier_type,'none')
            dsecond=sig.ax.*dx2-sig.aalpha.*sig.alpha1_dot+sig.Gamma.*sig.F_dot;
        else
            dsecond=dx2;
        end
        dy=[sc.g1.*x2+f1; dsecond; sig.ds1; sig.ds2; sig.dtheta1; sig.dtheta2];
        if strcmpi(p.alpha_dot_regressor_mode,'command_filter')
            dy=[dy;sig.dalpha_filter];
        end
        if p.state_derivative_filter
            dy=[dy;sig.dx2_filter];
        end
        if p.direction_estimator
            dy=[dy;sig.dx1_filter;sig.deta1;sig.deta2];
        end
        if p.second_layer_safety_flip
            dy=[dy;sig.dsafety_mode2];
        end
        if p.envelope_filter_hold
            dy=[dy;sig.dF_filter];
        end
        h=1e-6*(1+norm(y))/max(1,norm(dy));
        sp=local_control(tk+h,y+h*dy,p,N,M,beta,b,H,ref_map,sc);
        sm=local_control(tk-h,y-h*dy,p,N,M,beta,b,H,ref_map,sc);
        arr.alpha1_dot(k,:) = ((sp.alpha1-sm.alpha1)/(2*h))';
        arr.F_dot(k,:) = ((sp.F-sm.F)/(2*h))';
        arr.chi_dot_dir(k,:) = ((sp.chi-sm.chi)/(2*h))';
    end
    if strcmp(p.q_mode,'raw')
        residual1=H*f1-b*sc.yd1(tk);
    else
        residual1=beta.*f1-beta.*ref_map*sc.yd1(tk);
    end
    for f={'e','q','alpha1','F','chi','u','ratio','N1','N2','Np1','Np2', ...
    'R1','R2','ax','aalpha','Gamma','Areg','Freg','bar1','bar2','drive1','drive2','ds1','ds2', ...
            'dtheta1','dtheta2','alpha_filter','dalpha_filter','virtual_guard_blend', ...
            'safety_flip2','N2_effective', ...
            'safety_trigger2','dsafety_mode2', ...
            'N1_effective','direction_blend1','direction_blend2', ...
            'x2_filter','dx2_filter','ratio_trend_est', ...
            'x1_filter','dx1_filter','eta1','eta2','deta1','deta2', ...
            'lock_blend1','lock_blend2', ...
            'sigma1_eff','sigma2_eff','c2_eff'}
        arr.(f{1})(k,:)=sig.(f{1})';
    end
    arr.f1(k,:)=f1'; arr.f2(k,:)=f2'; arr.edot(k,:)=edot';
    arr.residual1(k,:)=residual1'; arr.yd(k)=sc.yd(tk); arr.ell(k)=sig.ell;
    if ~p.known_direction && ~strcmp(p.q_mode,'raw')
        lhs=sig.e'*H*diag(sc.g1)*sig.alpha1;

        rhs=sum(sc.g1.*beta.*sig.q.*sig.alpha1);
        arr.graph_identity_error(k)=lhs-rhs;
    end
end

if ~p.audit_detail
    arr.alpha1_dot=local_gradient(arr.alpha1,t);
    arr.F_dot=local_gradient(arr.F,t);
    arr.chi_dot_dir=local_gradient(arr.chi,t);
end
arr.chi_dot_fd=local_gradient(arr.chi,t);
arr.phi2=arr.ax.*arr.f2-arr.aalpha.*arr.alpha1_dot+arr.Gamma.*arr.F_dot;
arr.chi_dot_rhs=arr.ax.*arr.u.*sc.g2'+arr.phi2;
arr.barrier_derivative_error=arr.chi_dot_dir-arr.chi_dot_rhs;

theta1_req=max(abs(arr.residual1)./max(arr.R1,1e-12),[],1)+1e-10;
theta2_req=max(abs(arr.phi2)./max(arr.R2,1e-12),[],1)+1e-10;
arr.theta1_audit=repmat(theta1_req,n,1);
arr.theta2_audit=repmat(theta2_req,n,1);

t1=arr.theta1_audit-arr.theta1; t2=arr.theta2_audit-arr.theta2;
arr.V=0.5*sum(arr.e.^2,2)+0.5*sum(arr.chi.^2,2) ...
    +sum(t1.^2,2)/(2*p.gamma1)+sum(t2.^2,2)/(2*p.gamma2);
arr.V_dot=sum(arr.e.*arr.edot,2)+sum(arr.chi.*arr.chi_dot_rhs,2) ...
    -sum(t1.*arr.dtheta1,2)/p.gamma1-sum(t2.*arr.dtheta2,2)/p.gamma2;

eta=p.young_eta(:)';
cq=p.proof_mu*p.c1-(beta'.*p.g1_upper').^2./(2*eta);
arr.cq=repmat(cq,n,1); arr.cchi=p.proof_mu*arr.c2_eff-eta/2;
if p.known_direction
    arr.lyapunov_rhs=nan(n,1); arr.lyapunov_residual=nan(n,1);
else
    search=sum((arr.N1.*sc.g1'+p.proof_mu).*(arr.ds1./arr.sigma1_eff),2) ...
        +sum((arr.N2.*sc.g2'+p.proof_mu).*(arr.ds2./arr.sigma2_eff),2);
    leakage=-p.theta_leak1/(2*p.gamma1)*sum(t1.^2,2) ...
        -p.theta_leak2/(2*p.gamma2)*sum(t2.^2,2);
    floor_term=arr.ell*(sum(theta1_req)+sum(theta2_req)) ...
        +p.theta_leak1/(2*p.gamma1)*sum(theta1_req.^2) ...
        +p.theta_leak2/(2*p.gamma2)*sum(theta2_req.^2);
    arr.lyapunov_rhs=search-sum(arr.cq.*arr.q.^2,2) ...
        -sum(arr.cchi.*arr.chi.^2,2)+leakage+floor_term;
    arr.lyapunov_residual=arr.V_dot-arr.lyapunov_rhs;
end
arr.solver_step=[NaN;diff(t)];
arr.event_distance=1-max(abs(arr.ratio),[],2);
arr.cap_distance=p.exp_clip-max(local_nussbaum_exponent([arr.s1,arr.s2],p),[],2);
arr.beta=beta';
end

function out = local_metrics(arr,p,M,beta,sc,te,ie)
T=p.T; t=arr.t;
tmask=t>=T/2;
if sum(tmask)<2, tmask=false(size(t)); tmask(max(1,end-10):end)=true; end
e_norm=vecnorm(arr.e,2,2); chi_norm=vecnorm(arr.chi,2,2);
max_ratio=max(abs(arr.ratio),[],'all');
reason='none';
if ~isempty(ie)
    reasons={'envelope_event','nussbaum_cap_event','state_guard','command_guard','ratio_guard'};
    reason=reasons{ie(end)};
elseif t(end)<T-1.5*max(p.max_step,p.dt)
    reason='solver_termination';
end
svals=svd(M);
out.success=strcmp(reason,'none') && t(end)>=T-1.5*p.max_step && max_ratio<1;
out.tracking_success=out.success && sqrt(mean(e_norm(tmask).^2))<0.10;
out.safety_success=max_ratio<1 && ~strcmp(reason,'envelope_event');
out.reason=reason; out.T_final=t(end); out.fail_time=NaN;
if ~isempty(te), out.fail_time=te(end); end
out.final_norm_e=norm(arr.e(end,:));
out.rms_e_steady=sqrt(mean(e_norm(tmask).^2));
out.mean_e_steady=mean(e_norm(tmask));
out.rms_chi_steady=sqrt(mean(chi_norm(tmask).^2));
out.max_ratio=max_ratio;
out.min_safety_margin=min(arr.F-abs(arr.x2),[],'all');
out.max_abs_u=max(abs(arr.u),[],'all');
out.mean_u_steady=mean(max(abs(arr.u(tmask,:)),[],2));
out.lambda_beta=svals(end)^2; out.cond_M=svals(1)/svals(end);
out.beta=beta'; out.topology=p.topology; out.scenario=p.scenario;
out.solver=p.solver; out.rel_tol=p.rel_tol; out.abs_tol=p.abs_tol; out.max_step=p.max_step;
out.max_s=max(abs([arr.s1,arr.s2]),[],'all');
out.max_theta1=max(arr.theta1,[],'all'); out.max_theta2=max(arr.theta2,[],'all');
out.terminal_theta1=arr.theta1(end,:); out.terminal_theta2=arr.theta2(end,:);
out.cap_triggered=strcmp(reason,'nussbaum_cap_event');
out.envelope_event=strcmp(reason,'envelope_event');
out.state_guard_triggered=strcmp(reason,'state_guard');
out.command_guard_triggered=strcmp(reason,'command_guard');
out.ratio_guard_triggered=strcmp(reason,'ratio_guard');
out.graph_identity_error=max(abs(arr.graph_identity_error));
scale=max(1,max(abs(arr.chi_dot_rhs),[],'all'));
out.barrier_derivative_rel_error=max(abs(arr.barrier_derivative_error),[],'all')/scale;
out.theta1_audit=arr.theta1_audit(1,:); out.theta2_audit=arr.theta2_audit(1,:);
if all(isnan(arr.lyapunov_residual))
    out.lyapunov_residual_max=NaN;
else
    edge=max(2,round(0.02/max(median(diff(t)),eps)));
    idx=(1+edge):max(1,numel(t)-edge);
    if isempty(idx), idx=1:numel(t); end
    out.lyapunov_residual_max=max(arr.lyapunov_residual(idx));
end
out.theorem_admissible=p.lam_floor>0 && p.theta_leak1>0 && p.theta_leak2>0 ...
    && ~p.known_direction && isinf(p.command_limit) && out.success ...
    && strcmp(p.barrier_type,'artanh') && ~strcmp(p.q_mode,'raw') ...
    && (~p.direction_estimator || (p.state_derivative_filter ...
        && strcmpi(p.direction_switch_mode,'smooth'))) ...
    && ~p.virtual_command_guard ...
    && ~p.second_layer_safety_flip;
out.g1=sc.g1'; out.g2=sc.g2';
end

function D = local_gradient(X,t)
D=zeros(size(X));
if numel(t)<2, return; end
for j=1:size(X,2), D(:,j)=gradient(X(:,j),t); end
end

function [t,Y,te,ye,ie] = local_rk4(rhs,event_fun,y0,T,dt,save_stride)
n=round(T/dt); save_every=max(1,round(save_stride/dt));
n_save=floor(n/save_every)+2; t=zeros(n_save,1); Y=zeros(n_save,numel(y0));
t(1)=0; Y(1,:)=y0'; y=y0; tk=0; ks=1;
te=[]; ye=[]; ie=[];
for k=1:n
    k1=rhs(tk,y); k2=rhs(tk+dt/2,y+dt*k1/2);
    k3=rhs(tk+dt/2,y+dt*k2/2); k4=rhs(tk+dt,y+dt*k3);
    yn=y+dt*(k1+2*k2+2*k3+k4)/6;
    tnew=tk+dt;
    [v,~,~]=event_fun(tnew,yn);
    hit=find(v<=0,1);
    if ~isempty(hit)
        te=tnew; ye=yn'; ie=hit;
        if abs(t(ks)-tnew)>10*eps(max(1,tnew)), ks=ks+1; t(ks)=tnew; Y(ks,:)=yn'; end
        t=t(1:ks); Y=Y(1:ks,:); return;
    end
    y=yn; tk=tnew;
    if mod(k,save_every)==0 || k==n, ks=ks+1; t(ks)=tk; Y(ks,:)=y'; end
end
t=t(1:ks); Y=Y(1:ks,:);
end

function [A,b,H,hdiag] = local_graph_mats(topology)
switch lower(topology)
    case 'chain', A=[0 1 0;1 0 1;0 1 0]; b=[1;0;0];
    case 'star', A=[0 1 1;1 0 0;1 0 0]; b=[1;0;0];
    case 'complete', A=ones(3)-eye(3); b=[1;0;0];
    case 'directed', A=[0 0 1;1 0 0;0 1 0]; b=[1;0;0];
    otherwise, error('Unknown topology: %s',topology);
end
H=diag(sum(A,2)+b)-A; hdiag=diag(H);
end

function sc = local_scenario(name,p)
sc.g1=[-1.2;1.2;-0.8]; sc.g2=[2.0;-2.0;1.5];
switch lower(name)
    case {'offshore_nominal','offshore_ss2','nominal'}
        sc.x10=[0.35;-0.15;0.25]; sc.x20=[0.82;-0.78;0.68]; amp=0.15; freq=0.4; nf=0.03;
    case {'offshore_rough','offshore_ss4','rough'}
        sc.x10=[0.40;-0.22;0.30]; sc.x20=[0.55;-0.45;0.42]; amp=0.20; freq=0.5; nf=0.06;
    case 'offshore_ss1'
        sc.x10=[0.30;-0.10;0.20]; sc.x20=[0.70;-0.65;0.55]; amp=0.10; freq=0.3; nf=0.015;
    case 'offshore_ss3'
        sc.x10=[0.40;-0.20;0.30]; sc.x20=[0.65;-0.55;0.50]; amp=0.18; freq=0.45; nf=0.05;
    case {'decaying_diagnostic','asymptotic_regulation','asymptotic'}
        sc.x10=[0.35;-0.15;0.25]; sc.x20=zeros(3,1); amp=0; freq=0; nf=0;
    case 'physical_yaw_closedloop'
        sc.x10=zeros(3,1); sc.x20=zeros(3,1); amp=0; freq=0; nf=0;
    otherwise, error('Unknown scenario: %s',name);
end

load_scale=1;
if isfield(p,'mc_seed') && p.mc_seed>0
    rng(p.mc_seed); mode=p.mc_mode;
    switch lower(mode)
        case 'in_domain'
            sc.g1=sc.g1.*(0.95+0.10*rand(3,1)); sc.g2=sc.g2.*(0.95+0.10*rand(3,1));
            sc.x10=sc.x10+0.03*(2*rand(3,1)-1); sc.x20=sc.x20+0.03*(2*rand(3,1)-1);
        case 'sign_only'
            sc.g1=abs(sc.g1).*sign(randn(3,1)); sc.g2=abs(sc.g2).*sign(randn(3,1));
        case 'gain_only'
            sc.g1=sc.g1.*(0.90+0.20*rand(3,1)); sc.g2=sc.g2.*(0.90+0.20*rand(3,1));
        case 'initial_only'
            sc.x10=sc.x10+0.10*(2*rand(3,1)-1); sc.x20=sc.x20+0.10*(2*rand(3,1)-1);
        case 'load_only', load_scale=1+0.50*rand();
        case {'combined_stress','stress'}
            sc.g1=abs(sc.g1).*sign(randn(3,1)).*(0.90+0.20*rand(3,1));
            sc.g2=abs(sc.g2).*sign(randn(3,1)).*(0.90+0.20*rand(3,1));
            sc.x10=sc.x10+0.10*(2*rand(3,1)-1); sc.x20=sc.x20+0.10*(2*rand(3,1)-1);
        otherwise, error('Unknown Monte-Carlo mode: %s',mode);
    end
end
if isfield(p,'sign_pattern') && ~isempty(p.sign_pattern)
    sp=p.sign_pattern(:); sc.g1=abs(sc.g1).*sp(1:3); sc.g2=abs(sc.g2).*sp(4:6);
end
if isfield(p,'g1_override'), sc.g1=p.g1_override(:); end
if isfield(p,'g2_override'), sc.g2=p.g2_override(:); end
if isfield(p,'x10_override'), sc.x10=p.x10_override(:); end
if isfield(p,'x20_override'), sc.x20=p.x20_override(:); end
if isfield(p,'load_scale_override'), load_scale=p.load_scale_override; end

sc.yd=@(t) amp*sin(freq*t);
sc.yd1=@(t) amp*freq*cos(freq*t);
sc.yd2=@(t) -amp*freq^2*sin(freq*t);
sc.yd3=@(t) -amp*freq^3*cos(freq*t);
if strcmpi(name,'physical_yaw_closedloop')
    sc.f1=@(x1,~,~) zeros(size(x1));
    sc.f2=@(x1,~,~) zeros(size(x1));
elseif any(strcmpi(name,{'decaying_diagnostic','asymptotic_regulation','asymptotic'}))
    decay=@(t) exp(-0.2*t);
    sc.f1=@(x1,~,t) decay(t).*(0.03*x1.^2+0.01*sin(0.3*t+(1:3)'));
    sc.f2=@(x1,x2,t) decay(t).*(0.03*sin(x1.*x2)+0.01*cos(0.2*t+(1:3)'));
else
    sc.f1=@(x1,~,t) load_scale*(nf*x1.^2+nf*x1+0.01*sin(0.3*t+(1:3)'));
    sc.f2=@(x1,x2,t) load_scale*(nf*sin(x1.*x2)+0.01*cos(0.2*t+(1:3)'));
end
end

function [N,Np] = local_nussbaum(s,p,layer)
if layer==1, eps_=p.eps1; else, eps_=p.eps2; end
idx=p.nussbaum_index(layer,:)'; k=p.nussbaum_base.^(idx-1);
if strcmpi(p.nussbaum_family,'yu2025')
    a=p.yu_a; b0=p.yu_b; c=p.yu_c; w=p.yu_omega;
    if ~(isscalar(a) && a>0 && isscalar(b0) && b0>=0 ...
            && isscalar(c) && c>0.5 && c<1 && isscalar(w) && w>0)
        error('yu2025 requires a>0, b>=0, 0.5<c<1, and omega>0.');
    end
    f=a*(s.^2+b0).^c;
    E=exp(min(f,p.exp_clip));
    fp=2*a*c*s.*(s.^2+b0).^(c-1);
    active=f<p.exp_clip-10*eps(p.exp_clip);
    if strcmpi(p.yu_phase,'sin')
        N=eps_*E.*sin(w*s);
        Np=eps_*E.*w.*cos(w*s);
        Np(active)=Np(active)+eps_*E(active).*fp(active).*sin(w*s(active));
    else
        N=eps_*E.*cos(w*s);
        Np=-eps_*E.*w.*sin(w*s);
        Np(active)=Np(active)+eps_*E(active).*fp(active).*cos(w*s(active));
    end
elseif strcmpi(p.nussbaum_family,'ma2025')
    [N,Np]=local_ma_nussbaum(s,idx,p,eps_,false);
elseif strcmpi(p.nussbaum_family,'ma2025_smooth')
    [N,Np]=local_ma_nussbaum(s,idx,p,eps_,true);
elseif strcmpi(p.nussbaum_family,'qiao2022')
    a=p.qiao_a(idx); b=p.qiao_b(idx); T=p.qiao_T(idx);
    sa=abs(s);
    level=max(1,ceil(log(1+(b-1).*sa./(pi*T))./log(b)-1e-12));
    stretch=b.^(level-1);
    left=pi*T.*(stretch-1)./(b-1);
    phase=(sa-left)./(T.*stretch);
    parity=(-1).^(level-1);
    N=eps_*a.*parity.*cos(phase);
    Np=-eps_*a.*parity.*sin(phase)./(T.*stretch).*sign(s);
elseif strcmpi(p.nussbaum_family,'zhou2026')

    alpha=p.nussbaum_alpha; beta=p.nussbaum_beta; gamma=p.nussbaum_gamma;
    scale=alpha.^gamma.*beta.^idx;
    w=1./scale;
    C=sqrt(alpha^2*(gamma+1).*beta.^(2*idx)+1)./scale;
    exponent=min(alpha*abs(s),p.exp_clip); E=exp(exponent);
    N=eps_*C.*E.*sin(w.*s);
    active=alpha*abs(s)<p.exp_clip-10*eps(p.exp_clip);
    Np=eps_*C.*E.*w.*cos(w.*s);
    Np(active)=Np(active)+eps_*C(active).*E(active).*alpha.*sign(s(active)).*sin(w(active).*s(active));
elseif strcmpi(p.nussbaum_family,'huang2018')
    alpha=p.nussbaum_alpha; beta=p.nussbaum_beta;
    w=beta.^(-idx); C=sqrt(alpha^2*beta.^(2*idx)+1)./beta.^idx;
    if isfield(p,'nussbaum_normalize_amplitude') && p.nussbaum_normalize_amplitude

        C=ones(size(C));
    end
    exponent=min(alpha*abs(s),p.exp_clip); E=exp(exponent);
    N=eps_*C.*E.*sin(w.*s);
    active=alpha*abs(s)<p.exp_clip-10*eps(p.exp_clip);
    Np=eps_*C.*E.*w.*cos(w.*s);
    Np(active)=Np(active)+eps_*C(active).*E(active).*alpha.*sign(s(active)).*sin(w(active).*s(active));
else
    s2=min(s.^2,p.exp_clip); sh=sinh(s2); ch=cosh(s2);
    N=eps_*(2*s.*sh.*cos(k.*s)-k.*(ch-1).*sin(k.*s));
    inside=s.^2<p.exp_clip-10*eps(p.exp_clip);
    Np=eps_*(2*sh.*cos(k.*s)-2*k.*s.*sh.*sin(k.*s)-k.^2.*(ch-1).*cos(k.*s));
    Np(inside)=eps_*(2*sh(inside).*cos(k(inside).*s(inside)) ...
        +4*s(inside).^2.*ch(inside).*cos(k(inside).*s(inside)) ...
        -4*k(inside).*s(inside).*sh(inside).*sin(k(inside).*s(inside)) ...
        -k(inside).^2.*(ch(inside)-1).*cos(k(inside).*s(inside)));
end
end

function exponent = local_nussbaum_exponent(s,p)
if strcmpi(p.nussbaum_family,'yu2025')
    exponent=p.yu_a*(s.^2+p.yu_b).^p.yu_c;
elseif any(strcmpi(p.nussbaum_family,{'ma2025','ma2025_smooth','qiao2022'}))
    exponent=zeros(size(s));
elseif any(strcmpi(p.nussbaum_family,{'huang2018','zhou2026'}))
    exponent=p.nussbaum_alpha*abs(s);
else
    exponent=s.^2;
end
end

function [scale,scale_dot]=local_deadzone_schedule(t,p)
tau=p.search_deadzone_ramp_time;
if ~(isscalar(tau) && tau>=0)
    error('search_deadzone_ramp_time must be a nonnegative scalar.');
end
if tau==0
    scale=1; scale_dot=0;
else
    r=max(t,0)/tau;
    scale=1-exp(-r^2);
    scale_dot=2*max(t,0)/tau^2*exp(-r^2);
end
end

function [blend,time_dot,N_sensitivity,error_sensitivity]= ...
        local_search_lock_blend(t,p,N,error_signal,layer)
time_field=sprintf('search_lock_time%d',layer);
transition_field=sprintf('search_lock_transition%d',layer);
lock_time=p.search_lock_time;
transition=p.search_lock_transition;
if isfield(p,time_field) && ~isnan(p.(time_field)), lock_time=p.(time_field); end
if isfield(p,transition_field) && ~isnan(p.(transition_field))
    transition=p.(transition_field);
end
if isinf(lock_time)
    blend=zeros(size(N)); time_dot=zeros(size(N)); N_sensitivity=zeros(size(N));
    error_sensitivity=zeros(size(N)); return;
end
if ~(isscalar(lock_time) && lock_time>=0 && isscalar(transition) && transition>0)
    error('Finite search lock requires nonnegative lock time and positive transition.');
end
r=(t-lock_time)/transition;
if r<=0
    time_blend=0; time_blend_dot=0;
elseif r>=1
    time_blend=1; time_blend_dot=0;
else
    time_blend=3*r^2-2*r^3;
    time_blend_dot=6*r*(1-r)/transition;
end
gain_min=p.search_lock_gain_min; gain_full=p.search_lock_gain_full;
gain_min_field=sprintf('search_lock_gain_min%d',layer);
gain_full_field=sprintf('search_lock_gain_full%d',layer);
if isfield(p,gain_min_field) && ~isnan(p.(gain_min_field))
    gain_min=p.(gain_min_field);
end
if isfield(p,gain_full_field) && ~isnan(p.(gain_full_field))
    gain_full=p.(gain_full_field);
end
if ~(isscalar(gain_min) && gain_min>=0 && isscalar(gain_full) && gain_full>gain_min)
    error('Search-lock gain thresholds must satisfy 0 <= min < full.');
end
gain=abs(N); confidence=zeros(size(N)); dc_dN=zeros(size(N));
mid=gain>gain_min & gain<gain_full; high=gain>=gain_full;
confidence(high)=1;
z=(gain(mid)-gain_min)/(gain_full-gain_min);
confidence(mid)=3*z.^2-2*z.^3;
dc_dgain=6*z.*(1-z)/(gain_full-gain_min);
dc_dN(mid)=dc_dgain.*sign(N(mid));

low=p.(sprintf('search_lock_error_low%d',layer));
high=p.(sprintf('search_lock_error_high%d',layer));
if isinf(low) && isinf(high)
    error_confidence=ones(size(N)); dc_derr=zeros(size(N));
else
    if ~(isscalar(low) && low>=0 && isscalar(high) && high>low)
        error('Search-lock error thresholds must satisfy 0 <= low < high.');
    end
    magnitude=abs(error_signal);
    error_confidence=ones(size(N)); dc_derr=zeros(size(N));
    mid_error=magnitude>low & magnitude<high;
    high_error=magnitude>=high;
    error_confidence(high_error)=0;
    ze=(magnitude(mid_error)-low)/(high-low);
    error_confidence(mid_error)=1-(3*ze.^2-2*ze.^3);
    dc_dmag=-6*ze.*(1-ze)/(high-low);
    dc_derr(mid_error)=dc_dmag.*sign(error_signal(mid_error));
end

monitor_low=p.search_lock_monitor_gain_low;
monitor_high=p.search_lock_monitor_gain_high;
monitor_low_field=sprintf('search_lock_monitor_gain_low%d',layer);
monitor_high_field=sprintf('search_lock_monitor_gain_high%d',layer);
if isfield(p,monitor_low_field) && ~isnan(p.(monitor_low_field))
    monitor_low=p.(monitor_low_field);
end
if isfield(p,monitor_high_field) && ~isnan(p.(monitor_high_field))
    monitor_high=p.(monitor_high_field);
end
if ~(isscalar(monitor_low) && monitor_low>=gain_full && ...
        isscalar(monitor_high) && monitor_high>monitor_low)
    error(['Search-lock monitor gain thresholds must satisfy gain_full <= ' ...
        'monitor_low < monitor_high.']);
end
reliable=zeros(size(N)); dr_dN=zeros(size(N));
mid_reliable=gain>monitor_low & gain<monitor_high;
high_reliable=gain>=monitor_high;
reliable(high_reliable)=1;
zr=(gain(mid_reliable)-monitor_low)/(monitor_high-monitor_low);
reliable(mid_reliable)=3*zr.^2-2*zr.^3;
dr_dgain=6*zr.*(1-zr)/(monitor_high-monitor_low);
dr_dN(mid_reliable)=dr_dgain.*sign(N(mid_reliable));

monitored_confidence=error_confidence+reliable.*(1-error_confidence);
dmonitored_dN=dr_dN.*(1-error_confidence);
dmonitored_derr=(1-reliable).*dc_derr;
blend=time_blend.*confidence.*monitored_confidence;
time_dot=time_blend_dot.*confidence.*monitored_confidence;
N_sensitivity=time_blend.*(dc_dN.*monitored_confidence ...
    +confidence.*dmonitored_dN);
error_sensitivity=time_blend.*confidence.*dmonitored_derr;
end

function sigma=local_search_sigma(t,p,layer,N)
base=p.(sprintf('sigma%d',layer));
sigma=base.*ones(N,1);
if layer~=1 || isnan(p.search_sigma1_late) || isinf(p.search_sigma1_time)
    return;
end
transition=p.search_sigma1_transition;
if ~(isscalar(transition) && transition>0)
    error('search_sigma1_transition must be positive.');
end
r=(t-p.search_sigma1_time)/transition;
if r<=0
    blend=0;
elseif r>=1
    blend=1;
else
    blend=3*r^2-2*r^3;
end
sigma=(base+(p.search_sigma1_late-base)*blend).*ones(N,1);
end

function c2=local_c2_schedule(t,p,N)
c2=p.c2.*ones(N,1);
if isnan(p.search_c2_late) || isinf(p.search_c2_time), return; end
transition=p.search_c2_transition;
if ~(isscalar(transition) && transition>0)
    error('search_c2_transition must be positive.');
end
r=(t-p.search_c2_time)/transition;
if r<=0
    blend=0;
elseif r>=1
    blend=1;
else
    blend=3*r^2-2*r^3;
end
c2=(p.c2+(p.search_c2_late-p.c2)*blend).*ones(N,1);
end

function [N,Np] = local_ma_nussbaum(s,idx,p,eps_,smooth_variant)

alpha=p.ma_alpha; beta=p.ma_beta; M=p.ma_amp;
if ~(isscalar(alpha) && alpha>1 && isscalar(beta) && beta>0 ...
        && isscalar(M) && M>0)
    error('ma2025 requires scalar ma_alpha>1, ma_beta>0, and ma_amp>0.');
end
N=zeros(size(s)); Np=zeros(size(s));
for ii=1:numel(s)
    x=abs(s(ii)); block=2^(idx(ii)-1); n=0; left=0;
    radius=local_ma_radius(alpha,beta,block,n);
    right=left+2*radius;
    while x>=right && isfinite(right)
        left=right; n=n+1;
        if n>10000
            error('ma2025 interval search exceeded 10000 blocks.');
        end
        radius=local_ma_radius(alpha,beta,block,n);
        right=left+2*radius;
    end
    tau=x-left;
    if isinf(radius)
        inner=tau;
    else
        inner=radius-abs(tau-radius);
    end
    parity=(-1)^n;
    if smooth_variant
        if isinf(radius)
            phase=0;
        else
            phase=pi*tau/(2*radius);
        end
        N(ii)=eps_*M*parity*sin(phase)^2;
        Np(ii)=eps_*M*parity*pi/(2*radius)*sin(2*phase)*sign(s(ii));
        continue;
    end
    denom=tanh(radius);
    N(ii)=eps_*M*parity*tanh(inner)/denom;
    if isinf(radius) || tau<radius
        dinner=1;
    elseif tau>radius
        dinner=-1;
    else

        dinner=0;
    end
    Np(ii)=eps_*M*parity*(1-tanh(inner)^2)/denom*dinner*sign(s(ii));
end
end

function radius = local_ma_radius(alpha,beta,block,n)
j=(block*n):(block*(n+1)-1);
logr=log(beta)+0.5*j.*(j-1)*log(alpha);
if max(logr)>log(realmax)-log(max(block,2))
    radius=Inf;
else
    radius=sum(exp(logr));
end
end

function [F,Falpha,Ffilter,Fcmd] = local_envelope(t,alpha,envelope_filter,p)
amplitude=sqrt(alpha.^2+p.envelope_eps_abs^2);
A=amplitude/p.r;
root=sqrt((p.rho-A).^2+p.envelope_eps_max^2);
Fcmd=0.5*(p.rho+A+root);
dFdA=0.5*(1+(A-p.rho)./root);
Fcmd_alpha=dFdA.*alpha./(p.r*amplitude);
if p.envelope_filter_hold
    held_value=min(envelope_filter,Fcmd+p.envelope_hold_gap_max);
    use_hold=t<p.envelope_hold_end_time & held_value>Fcmd ...
        & envelope_filter>p.rho+p.envelope_hold_activation;
    cap_active=use_hold & envelope_filter>Fcmd+p.envelope_hold_gap_max;
    free_hold=use_hold & ~cap_active;
    F=Fcmd; F(use_hold)=held_value(use_hold);
    Falpha=Fcmd_alpha; Falpha(free_hold)=0;
    Ffilter=zeros(size(alpha)); Ffilter(free_hold)=1;
else
    F=Fcmd;
    Falpha=Fcmd_alpha;
    Ffilter=zeros(size(alpha));
end
end

function [B,dB] = local_barrier(z,p)
if strcmp(p.barrier_type,'none'), B=z; dB=ones(size(z)); return; end
zz=max(min(z./p.barrier_margin,p.ratio_clip),-p.ratio_clip);
B=atanh(zz); dB=1./(p.barrier_margin*(1-zz.^2));
end

function [w,dwdv] = local_unit_smooth(v,low,high)
x=max(0,min(1,(v-low)./(high-low)));
w=x.^2.*(3-2*x);
dwdv=zeros(size(v));
inside=v>low & v<high;
dwdv(inside)=6*x(inside).*(1-x(inside))./(high-low);
end

function [Neff,blend,blend_eta] = local_direction_override(N,eta,p,layer)
if ~p.direction_estimator
    Neff=N; blend=zeros(size(N)); blend_eta=zeros(size(N)); return;
end
high=p.(sprintf('direction_confidence%d_full',layer));
if strcmpi(p.direction_switch_mode,'relay')
    blend=double(abs(eta)>=high);
    blend_eta=zeros(size(eta));
else
    low=p.(sprintf('direction_confidence%d_min',layer));
    [blend,dwdabs]=local_unit_smooth(abs(eta),low,high);
    blend_eta=dwdabs.*sign(eta);
end
identified=-sign(eta); identified(identified==0)=1;
Neff=(1-blend).*N+blend.*identified;
end

function z = local_barrier_inverse(B,p)
if strcmp(p.barrier_type,'none'), z=B; return; end
z=p.barrier_margin*tanh(B);
end

function [d,dp] = local_damp(z,p,layer)
if ~p.nonlinear_damping, d=zeros(size(z)); dp=zeros(size(z)); return; end
if layer==1, ka=p.kfa1; kb=p.kfb1; kp=p.probe_gain1;
else, ka=p.kfa2; kb=p.kfb2; kp=p.probe_gain2; end
S=z.^2+p.eps_q^2; m=(1-p.qa)/2; n=(p.qb-1)/2;
d=ka*z./S.^m+kb*S.^n.*z+kp*tanh(z/p.probe_eps);
dp=ka*(S.^(-m)-2*m*z.^2.*S.^(-m-1)) ...
    +kb*(S.^n+2*n*z.^2.*S.^(n-1)) ...
    +kp/p.probe_eps*(1-tanh(z/p.probe_eps).^2);
end

function beta = local_opt_beta(H,hdiag,p)
N=size(H,1); if N~=3, beta=hdiag(:); return; end
grid=p.beta_min:p.beta_step:p.beta_max; best_obj=-Inf; beta=hdiag(:);
for b1=grid
    for b2=grid
        b3=N*p.beta_bar-b1-b2;
        if b3<p.beta_min-1e-12 || b3>p.beta_max+1e-12, continue; end
        bt=[b1;b2;b3]; lam=min(eig(H*diag(1./bt.^2)*H'));
        obj=lam-p.beta_rho_u*sum((bt-hdiag(:)).^2);
        if obj>best_obj, best_obj=obj; beta=bt; end
    end
end
end
