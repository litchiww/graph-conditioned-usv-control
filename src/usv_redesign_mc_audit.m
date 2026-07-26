function T=usv_redesign_mc_audit(root_dir,n_runs)

if nargin<1, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
if nargin<2, n_runs=30; end
cfg=usv_config();
schema_version=30;
basep=usv_redesign_candidate(cfg.get_params('C2_full'));
signature=jsonencode(basep);
checkpoint=fullfile(root_dir,'data','sign_redesign_mc_checkpoint.mat');
rows=[];
if exist(checkpoint,'file')
    S=load(checkpoint);
    if isfield(S,'schema_version') && S.schema_version==schema_version && ...
            isfield(S,'signature') && strcmp(S.signature,signature) && ...
            isfield(S,'n_runs') && S.n_runs==n_runs && isfield(S,'rows')
        rows=S.rows;
    end
end
start_id=height(rows)+1;
for run_id=start_id:n_runs
    p=basep;
    p.mc_seed=7000+run_id; p.mc_mode='combined_stress';
    p.audit_detail=false; p.save_stride=0.02; p.max_wall_time=300;
    timer=tic; [o,a]=usv_simulate(p); runtime=toc(timer);
    truth=sign([o.g1 o.g2])';
    estimate=[sign(a.eta1(end,:)) sign(a.eta2(end,:))]';
    sign_match=all(estimate==truth);
    accepted=o.tracking_success && o.success && o.max_ratio<p.max_ratio_guard ...
        && sign_match;
    row=table(run_id,p.mc_seed,string(sprintf('%+d',truth)), ...
        string(sprintf('%+d',estimate)),sign_match,o.success,string(o.reason), ...
        o.T_final,o.rms_e_steady,o.max_ratio,o.max_abs_u,o.max_s, ...
        min(abs([a.eta1(end,:) a.eta2(end,:)])), ...
        min([a.direction_blend1(end,:) a.direction_blend2(end,:)]), ...
        accepted,runtime, ...
        'VariableNames',{'run_id','seed','true_signs','estimated_signs', ...
        'sign_match','completed','reason','T_final','rms_e','max_ratio', ...
        'peak_u','max_s','min_abs_eta','min_confidence_blend','accepted','runtime_s'});
    rows=[rows;row];
    save(checkpoint,'rows','schema_version','signature','n_runs','-v7.3');
    writetable(rows,fullfile(root_dir,'data','sign_redesign_mc_partial.csv'));
    fprintf('MC %02d/%02d %s RMS %.4f ratio %.4f peak %.1f signs %d\n', ...
        run_id,n_runs,o.reason,o.rms_e_steady,o.max_ratio,o.max_abs_u,sign_match);
    drawnow;
end
T=rows;
writetable(T,fullfile(root_dir,'data','sign_redesign_mc.csv'));
save(fullfile(root_dir,'data','sign_redesign_mc.mat'), ...
    'rows','schema_version','signature','n_runs','-v7.3');
end
