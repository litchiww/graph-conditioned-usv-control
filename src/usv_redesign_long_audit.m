function T=usv_redesign_long_audit(root_dir)

if nargin<1, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
cfg=usv_config();
basep=usv_redesign_candidate(cfg.get_params('C2_full'));
schema_version=30;
worst_code=3;
sign_file=fullfile(root_dir,'data','sign_redesign_patterns.csv');
if exist(sign_file,'file')
    sign_rows=readtable(sign_file);
    sign_rows=sign_rows(logical(sign_rows.accepted),:);
    [~,idx]=max(sign_rows.max_ratio);
    worst_code=sign_rows.pattern_id(idx);
end
specs={ ...
    'nominal',      42, 'chain',   '',        1.0; ...
    'worst_ratio', worst_code, 'chain', '',        1.0; ...
    'star',         42, 'star',    '',        1.0; ...
    'initial_edge', 42, 'chain',   'initial', 1.0; ...
    'load_1p5',     42, 'chain',   '',        1.5};
signature=jsonencode(struct('basep',basep,'specs',{specs}));
checkpoint=fullfile(root_dir,'data','sign_redesign_long_checkpoint.mat');
rows=table();
if exist(checkpoint,'file')
    S=load(checkpoint);
    if isfield(S,'schema_version') && S.schema_version==schema_version && ...
            isfield(S,'signature') && strcmp(S.signature,signature) && isfield(S,'rows')
        rows=S.rows;
    end
end

pending=true(size(specs,1),1);
if ~isempty(rows)
    for k=1:size(specs,1)
        pending(k)=~any(rows.case_name==string(specs{k,1}));
    end
end
indices=find(pending);
if ~isempty(indices) && license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if isempty(pool), pool=parpool('Processes',min(3,numel(indices))); end
    futures(1,numel(indices))=parallel.FevalFuture;
    for k=1:numel(indices)
        futures(k)=parfeval(pool,@local_run_long,1,basep,specs(indices(k),:));
    end
    for k=1:numel(futures)
        [~,row]=fetchNext(futures);
        rows=[rows;row];
        save(checkpoint,'rows','schema_version','signature','-v7.3');
        writetable(rows,fullfile(root_dir,'data','sign_redesign_long_partial.csv'));
        local_report(row);
    end
else
    for k=indices'
        row=local_run_long(basep,specs(k,:));
        rows=[rows;row];
        save(checkpoint,'rows','schema_version','signature','-v7.3');
        writetable(rows,fullfile(root_dir,'data','sign_redesign_long_partial.csv'));
        local_report(row);
    end
end
T=sortrows(rows,'case_name');
writetable(T,fullfile(root_dir,'data','sign_redesign_long.csv'));
save(fullfile(root_dir,'data','sign_redesign_long.mat'), ...
    'rows','schema_version','signature','-v7.3');
end

function row=local_run_long(basep,spec)
p=basep;
code=spec{2}; p.sign_pattern=2*double(bitget(uint8(code),1:6))'-1;
p.topology=spec{3};
if strcmp(spec{4},'initial')
    p.x10_override=[0.45;-0.25;0.35];
    p.x20_override=[0.85;-0.81;0.72];
end
p.load_scale_override=spec{5};
p.T=600; p.audit_detail=false; p.save_stride=0.05; p.max_wall_time=1800;
timer=tic; [o,a]=usv_simulate(p); runtime=toc(timer);
truth=sign([o.g1 o.g2])';
estimate=[sign(a.eta1(end,:)) sign(a.eta2(end,:))]';
sign_match=all(estimate==truth);
accepted=o.tracking_success && o.success && o.max_ratio<p.max_ratio_guard ...
    && sign_match;
row=table(string(spec{1}),code,string(p.topology),o.success, ...
    string(o.reason),o.T_final,o.rms_e_steady,o.max_ratio,o.max_abs_u, ...
    o.max_s,sign_match,accepted,runtime, ...
    'VariableNames',{'case_name','pattern_id','topology','completed', ...
    'reason','T_final','rms_e','max_ratio','peak_u','max_s', ...
    'sign_match','accepted','runtime_s'});
end

function local_report(row)
fprintf('LONG %s %s RMS %.4f ratio %.4f peak %.1f signs %d\n', ...
    row.case_name,row.reason,row.rms_e,row.max_ratio,row.peak_u,row.sign_match);
drawnow;
end
