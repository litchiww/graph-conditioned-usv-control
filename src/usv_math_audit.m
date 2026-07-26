function report = usv_math_audit(arr, p, out_file)

rng(2401);
s = -1.25+2.5*rand(200,1);
h = 1e-20;
primitive_err = 0;
for ell = 1:6
    eps_ = p.eps1;
    if ell > 3, eps_ = p.eps2; end
    if strcmpi(p.nussbaum_family,'huang2018')
        alpha=p.nussbaum_alpha; beta=p.nussbaum_beta; w=beta^(-ell);
        C=sqrt(alpha^2*beta^(2*ell)+1)/beta^ell;
        if isfield(p,'nussbaum_normalize_amplitude') && p.nussbaum_normalize_amplitude
            C=1;
        end

        hr=1e-6;
        primitive=@(x) eps_*C.*(exp(alpha*abs(x)).*(alpha*sin(w*abs(x))-w*cos(w*abs(x)))+w) ...
            /(alpha^2+w^2);
        dpsi=(primitive(s+hr)-primitive(s-hr))/(2*hr);
        nfun=eps_*C.*exp(alpha*abs(s)).*sin(w*s);
    else
        k = p.nussbaum_base^(ell-1);
        z = s+1i*h;
        psi = eps_*(cosh(z.^2)-1).*cos(k*z);
        dpsi = imag(psi)/h;
        nfun = eps_*(2*s.*sinh(s.^2).*cos(k*s) ...
            -k*(cosh(s.^2)-1).*sin(k*s));
    end
    primitive_err = max(primitive_err,max(abs(dpsi-nfun)./max(1,abs(nfun))));
end

theta1 = arr.theta1_audit(1,:);
theta2 = arr.theta2_audit(1,:);
reg1 = max(abs(arr.residual1)./max(arr.R1.*theta1,1e-14),[],'all');
reg2 = max(abs(arr.phi2)./max(arr.R2.*theta2,1e-14),[],'all');
graph_err = max(abs(arr.graph_identity_error));
barrier_scale = max(1,max(abs(arr.chi_dot_rhs),[],'all'));
barrier_err = max(abs(arr.barrier_derivative_error),[],'all')/barrier_scale;
lyap_max = max(arr.lyapunov_residual,[],'omitnan');

report = struct();
report.nussbaum_primitive_relerr = primitive_err;
report.graph_identity_abserr = graph_err;
report.barrier_derivative_relerr = barrier_err;
report.regressor1_max_ratio = reg1;
report.regressor2_max_ratio = reg2;
report.lyapunov_residual_max = lyap_max;
report.pass_primitive = primitive_err < 1e-9;
report.pass_graph = graph_err < 1e-10;
report.pass_barrier = barrier_err < 1e-7;
report.pass_regressor = reg1 <= 1+1e-9 && reg2 <= 1+1e-9;
report.pass_lyapunov = lyap_max <= 1e-6;
report.pass_all = report.pass_primitive && report.pass_graph && ...
    report.pass_barrier && report.pass_regressor && report.pass_lyapunov;

if nargin >= 3 && ~isempty(out_file)
    fid = fopen(out_file,'w');
    if fid < 0, error('Could not write audit file: %s',out_file); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'test,value,tolerance,pass\n');
    fprintf(fid,'nussbaum_primitive,%.12g,1e-9,%d\n',primitive_err,report.pass_primitive);
    fprintf(fid,'graph_cancellation,%.12g,1e-10,%d\n',graph_err,report.pass_graph);
    fprintf(fid,'barrier_derivative,%.12g,1e-7,%d\n',barrier_err,report.pass_barrier);
    fprintf(fid,'composite_regressor_layer1,%.12g,1,%d\n',reg1,reg1<=1+1e-9);
    fprintf(fid,'composite_regressor_layer2,%.12g,1,%d\n',reg2,reg2<=1+1e-9);
    fprintf(fid,'lyapunov_inequality,%.12g,1e-6,%d\n',lyap_max,report.pass_lyapunov);
end
end
