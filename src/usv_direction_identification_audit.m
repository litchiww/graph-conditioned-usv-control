function T=usv_direction_identification_audit(root_dir)

if nargin<1, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
cfg=usv_config(); rows=[];
for code=0:63
    p=usv_redesign_candidate(cfg.get_params('C2_full'));
    p.sign_pattern=2*double(bitget(uint8(code),1:6))'-1;
    p.T=0.50; p.audit_detail=false; p.save_stride=0.001;
    [o,a]=usv_simulate(p);
    eta=[a.eta1 a.eta2];
    blend=[a.direction_blend1 a.direction_blend2];
    estimate=sign(eta(end,:))';
    truth=p.sign_pattern(:);
    active=blend>1e-12;
    mismatch=active & sign(eta)~=repmat(truth',size(eta,1),1);
    active_time=NaN(1,numel(truth));
    settled_time=NaN(size(active_time));
    for j=1:numel(truth)
        idx=find(active(:,j),1,'first');
        if ~isempty(idx), active_time(j)=a.t(idx); end
        idx=find(blend(:,j)>=0.99,1,'first');
        if ~isempty(idx), settled_time(j)=a.t(idx); end
    end
    row=table(code,string(sprintf('%+d',truth)),string(sprintf('%+d',estimate)), ...
        all(estimate==truth),~any(mismatch,'all'),max(active_time,[],'omitnan'), ...
        max(settled_time,[],'omitnan'),min(abs(eta(end,:))),min(blend(end,:)), ...
        o.max_ratio,o.max_abs_u,string(o.reason), ...
        'VariableNames',{'pattern_id','true_signs','estimated_signs','all_match', ...
        'correct_whenever_active','latest_activation_s','latest_99pct_blend_s', ...
        'min_abs_eta','min_terminal_blend','max_ratio','peak_u','reason'});
    rows=[rows;row];
end
T=rows;
writetable(T,fullfile(root_dir,'data','sign_direction_identification_audit.csv'));
if ~all(T.all_match & T.correct_whenever_active)
    error('Direction-identification audit found at least one active sign mismatch.');
end
end
