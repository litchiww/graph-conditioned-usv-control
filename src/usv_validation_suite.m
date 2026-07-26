function result = usv_validation_suite(stage, root_dir)

if nargin < 2, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
data_dir=fullfile(root_dir,'data'); case_dir=fullfile(root_dir,'cases');
if ~exist(data_dir,'dir'), mkdir(data_dir); end
if ~exist(case_dir,'dir'), mkdir(case_dir); end
cfg=usv_config();

switch lower(stage)
    case 'math'
        p=cfg.get_params('C2_full'); p.audit_detail=true;
        [out,arr]=usv_simulate(p);
        report=usv_math_audit(arr,p,fullfile(data_dir,'math_audit.csv'));
        local_gain_certificate(out,p,fullfile(data_dir,'gain_certificate.csv'));
        save(fullfile(case_dir,'C2_full.mat'),'out','arr','p','-v7.3');
        if ~report.pass_all
            error('SubmissionGate:MathAuditFailed', ...
                'Mathematical audit failed; inspect data/math_audit.csv.');
        end
        result=report;

    case 'solver'
        names=cfg.case_list(); dts=[5e-3,2.5e-3,1.25e-3]; rows=[];
        for i=1:numel(names)
            p=cfg.get_params(names{i}); p.audit_detail=false; p.save_stride=0.05;
            pref=p; pref.solver='ode15s';
            [oref,~]=usv_simulate(pref);
            for d=dts
                pt=p; pt.solver='rk4'; pt.dt=d; pt.save_stride=0.05;
                [o,~]=usv_simulate(pt);
                row=struct('case_name',names{i},'solver','RK4','dt',d, ...
                    'rms_e',o.rms_e_steady,'max_ratio',o.max_ratio,'peak_u',o.max_abs_u, ...
                    'max_s',o.max_s,'event_same',strcmp(o.reason,oref.reason), ...
                    'rms_rel_pct',100*abs(o.rms_e_steady-oref.rms_e_steady)/max(oref.rms_e_steady,eps), ...
                    'ratio_abs_diff',abs(o.max_ratio-oref.max_ratio), ...
                    'peak_rel_pct',100*abs(o.max_abs_u-oref.max_abs_u)/max(oref.max_abs_u,eps), ...
                    'accepted',false);
                row.accepted=row.rms_rel_pct<2 && row.peak_rel_pct<2 && row.event_same;
                rows=[rows;row];
            end
        end
        T=struct2table(rows); writetable(T,fullfile(data_dir,'solver_audit.csv'));
        save(fullfile(data_dir,'solver_audit.mat'),'rows','-v7.3'); result=rows;

    case {'signs','signs_redesign','signs_redesign_precheck','signs_redesign_scan'}
        basep=cfg.get_params('C2_full');
        is_redesign=startsWith(lower(stage),'signs_redesign');
        precheck_only=strcmpi(stage,'signs_redesign_precheck');
        scan_only=strcmpi(stage,'signs_redesign_scan');
        if is_redesign
            basep=local_redesign_params(basep);
            file_tag='sign_redesign'; schema_version=30;
        else
            file_tag='sign'; schema_version=2;
        end
        signature=jsonencode(basep);
        if ~scan_only
            if is_redesign, probe_codes=[0,2,3,4,7,42,57,59,63];
            else, probe_codes=[0,42,57,59,63]; end
            probe_rows=[];
            for code=probe_codes
                signs=2*double(bitget(uint8(code),1:6))'-1;
                p=basep; p.sign_pattern=signs; p.audit_detail=false;
                if is_redesign, p.save_stride=0.01; else, p.save_stride=0.05; end
                if is_redesign, p.max_wall_time=240; else, p.max_wall_time=120; end
                timer=tic; [o,~]=usv_simulate(p);
                probe_rows=[probe_rows;struct('pattern_id',code,'signs',sprintf('%+d',signs), ...
                    'completed',o.success,'reason',o.reason,'T_final',o.T_final, ...
                    'rms_e',o.rms_e_steady,'max_ratio',o.max_ratio,'peak_u',o.max_abs_u, ...
                    'max_s',o.max_s,'tracking_success',o.tracking_success, ...
                    'margin_success',o.success && o.max_ratio<basep.max_ratio_guard, ...
                    'accepted',o.tracking_success && o.success && o.max_ratio<basep.max_ratio_guard, ...
                    'runtime_s',toc(timer))];
            end
            writetable(struct2table(probe_rows),fullfile(data_dir,[file_tag '_gate_precheck.csv']));
            if ~all([probe_rows.accepted])
                error(['Fixed-sign hard gate failed in the representative precheck. ' ...
                    'Acceptance requires completion, RMS < 0.10, and continuous ratio < 0.98. ' ...
                    'See the gate precheck CSV; stale 64-pattern outputs are not valid.']);
            end
            if precheck_only
                result=probe_rows;
                return;
            end
        end
        checkpoint=fullfile(data_dir,[file_tag '_patterns_checkpoint.mat']);
        rows=repmat(struct('pattern_id',0,'signs','','completed',false,'reason','pending', ...
            'T_final',NaN,'rms_e',NaN,'max_ratio',NaN,'peak_u',NaN,'max_s',NaN,'max_theta',NaN, ...
            'envelope_event',false,'cap_triggered',false,'state_guard',false, ...
            'command_guard',false,'ratio_guard',false,'tracking_success',false,'margin_success',false, ...
            'accepted',false,'theorem_admissible',false,'runtime_s',NaN),64,1);
        if exist(checkpoint,'file')
            S=load(checkpoint);
            if isfield(S,'schema_version') && S.schema_version==schema_version && ...
                    isfield(S,'signature') && strcmp(S.signature,signature)
                rows=S.rows;
            end
        end
        for code=0:63
            if ~strcmp(rows(code+1).reason,'pending'), continue; end
            signs=2*double(bitget(uint8(code),1:6))'-1;
            p=basep; p.sign_pattern=signs; p.audit_detail=false;
            if ~is_redesign
                p.max_abs_state=Inf; p.max_abs_u=Inf;
            end
            if is_redesign, p.save_stride=0.01; else, p.save_stride=0.05; end
            if is_redesign, p.max_wall_time=240; else, p.max_wall_time=120; end
            timer=tic;
            [o,~]=usv_simulate(p);
            k=code+1; rows(k).pattern_id=code; rows(k).signs=sprintf('%+d',signs);
            rows(k).completed=o.success; rows(k).reason=o.reason; rows(k).T_final=o.T_final;
            rows(k).rms_e=o.rms_e_steady;
            rows(k).max_ratio=o.max_ratio; rows(k).peak_u=o.max_abs_u; rows(k).max_s=o.max_s;
            rows(k).max_theta=max(o.max_theta1,o.max_theta2);
            rows(k).envelope_event=o.envelope_event; rows(k).cap_triggered=o.cap_triggered;
            rows(k).state_guard=o.state_guard_triggered;
            rows(k).command_guard=o.command_guard_triggered;
            rows(k).ratio_guard=o.ratio_guard_triggered;
            rows(k).tracking_success=o.tracking_success;
            rows(k).margin_success=o.success && o.max_ratio<basep.max_ratio_guard;
            rows(k).accepted=o.tracking_success && rows(k).margin_success;
            rows(k).theorem_admissible=o.theorem_admissible;
            rows(k).runtime_s=toc(timer);
            save(checkpoint,'rows','schema_version','signature','-v7.3');
            T=struct2table(rows(~strcmp({rows.reason},'pending')));
            writetable(T,fullfile(data_dir,[file_tag '_patterns_partial.csv']));
            fprintf('sign %02d/63 %s ratio %.4f RMS %.4f runtime %.1fs\n', ...
                code,rows(k).reason,o.max_ratio,o.rms_e_steady,rows(k).runtime_s); drawnow;
            if is_redesign && ~rows(k).accepted
                error(['Redesign scan stopped at the first failed pattern. ' ...
                    'The checkpoint and partial CSV contain the diagnostic row.']);
            end
        end
        T=struct2table(rows); writetable(T,fullfile(data_dir,[file_tag '_patterns.csv']));
        save(fullfile(data_dir,[file_tag '_patterns.mat']),'rows','schema_version','signature','-v7.3'); result=rows;

    case 'sensitivity'
        rows=[]; zetas=[0.001,0.0025,0.005,0.01,0.02];
        lambdas=[1e-4,5e-4,1e-3,5e-3];
        for z=zetas
            p=cfg.get_params('C2_full'); p.theta_leak1=z; p.theta_leak2=z;
            p.audit_detail=false; p.save_stride=0.05; [o,arr]=usv_simulate(p);
            rows=[rows;local_sensitivity_row('leakage',z,o,arr,p)];
        end
        for lam=lambdas
            p=cfg.get_params('C2_full'); p.lam_floor=lam;
            p.audit_detail=false; p.save_stride=0.05; [o,arr]=usv_simulate(p);
            rows=[rows;local_sensitivity_row('denominator',lam,o,arr,p)];
        end
        T=struct2table(rows); writetable(T,fullfile(data_dir,'parameter_sensitivity.csv'));
        save(fullfile(data_dir,'parameter_sensitivity.mat'),'rows','-v7.3'); result=rows;

    case 'long'
        specs={ ...
            'nominal','C2_full',[]; ...
            'star','C3_star',[]; ...
            'largest_initial','C2_full','initial'; ...
            'largest_load','C2_full','load'; ...
            'worst_sign','C2_full','sign'};
        worst_sign=[];
        sf=fullfile(data_dir,'sign_patterns.mat');
        if exist(sf,'file')
            S=load(sf); ok=[S.rows.completed]; vals=[S.rows.rms_e]; vals(~ok)=-Inf;
            [~,ix]=max(vals); code=S.rows(ix).pattern_id;
            worst_sign=2*double(bitget(uint8(code),1:6))'-1;
        end
        rows=[];
        for k=1:size(specs,1)
            p=cfg.get_params(specs{k,2}); p.T=600; p.audit_detail=false; p.save_stride=0.05;
            tag=specs{k,3};
            if strcmp(tag,'initial')
                p.x10_override=[0.45;-0.25;0.35]; p.x20_override=[0.85;-0.81;0.72];
            elseif strcmp(tag,'load'), p.load_scale_override=1.5;
            elseif strcmp(tag,'sign') && ~isempty(worst_sign), p.sign_pattern=worst_sign;
            end
            [o,arr]=usv_simulate(p);
            row=struct('case_name',specs{k,1},'completed',o.success,'reason',o.reason, ...
                'late_rms',o.rms_e_steady,'max_ratio',o.max_ratio,'peak_u',o.max_abs_u, ...
                'max_s',o.max_s,'terminal_s',max(abs([arr.s1(end,:),arr.s2(end,:)])), ...
                'max_theta',max([o.max_theta1,o.max_theta2]), ...
                'terminal_theta',max([o.terminal_theta1,o.terminal_theta2]), ...
                'cap_triggered',o.cap_triggered);
            rows=[rows;row];
            out=o;
            save(fullfile(case_dir,['long_' specs{k,1} '.mat']),'out','arr','p','-v7.3');
        end
        T=struct2table(rows); writetable(T,fullfile(data_dir,'long_horizon.csv'));
        save(fullfile(data_dir,'long_horizon.mat'),'rows','-v7.3'); result=rows;

    otherwise
        error('Unknown validation stage: %s',stage);
end
end

function p=local_redesign_params(p)

p=usv_redesign_candidate(p);
end

function row=local_sensitivity_row(kind,value,o,arr,p)
du=diff(arr.u)./max(diff(arr.t),eps);
theta_terminal=max([arr.theta1(end,:),arr.theta2(end,:)]);
theta_drift=max(abs([arr.theta1(end,:)-arr.theta1(round(end/2),:), ...
    arr.theta2(end,:)-arr.theta2(round(end/2),:)]));
nu=min([2*min(arr.cq(1,:))*o.lambda_beta,2*min(arr.cchi(1,:)), ...
    p.theta_leak1,p.theta_leak2]);
row=struct('parameter',kind,'value',value,'rms_e',o.rms_e_steady, ...
    'peak_u',o.max_abs_u,'control_smoothness',sqrt(mean(du.^2,'all')), ...
    'estimate_drift',theta_drift,'terminal_estimate',theta_terminal, ...
    'max_ratio',o.max_ratio,'nu',nu,'completed',o.success);
end

function local_gain_certificate(out,p,out_file)
beta=out.beta(:); eta=p.young_eta(:);
cq=p.proof_mu*p.c1-(beta.*p.g1_upper).^2./(2*eta);
cchi=p.proof_mu*p.c2-eta/2;
nu=min([2*min(cq)*out.lambda_beta,2*min(cchi),p.theta_leak1,p.theta_leak2]);
T=table((1:3)',beta,p.g1_upper,repmat(p.proof_mu,3,1),eta, ...
    repmat(p.theta_leak1,3,1), ...
    repmat(p.theta_leak2,3,1),cq,cchi,repmat(out.lambda_beta,3,1),repmat(nu,3,1), ...
    'VariableNames',{'agent','beta','g1_upper','mu','young_weight', ...
    'leakage1','leakage2','c_q','c_chi','lambda_beta','nu'});
writetable(T,out_file);
end
