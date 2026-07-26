function rows=usv_redesign_sign_parallel(root_dir,n_workers)

if nargin<1, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
if nargin<2, n_workers=4; end
data_dir=fullfile(root_dir,'data');
cfg=usv_config();
basep=usv_redesign_candidate(cfg.get_params('C2_full'));
schema_version=30;
signature=jsonencode(basep);
checkpoint=fullfile(data_dir,'sign_redesign_patterns_checkpoint.mat');
rows=repmat(struct('pattern_id',0,'signs','','completed',false,'reason','pending', ...
    'T_final',NaN,'rms_e',NaN,'max_ratio',NaN,'peak_u',NaN,'max_s',NaN,'max_theta',NaN, ...
    'envelope_event',false,'cap_triggered',false,'state_guard',false, ...
    'command_guard',false,'ratio_guard',false,'tracking_success',false,'margin_success',false, ...
    'accepted',false,'theorem_admissible',false,'runtime_s',NaN),64,1);
if exist(checkpoint,'file')
    S=load(checkpoint);
    if isfield(S,'schema_version') && S.schema_version==schema_version && ...
            isfield(S,'signature') && strcmp(S.signature,signature) && isfield(S,'rows')
        rows=S.rows;
    end
end

pending=find(strcmp({rows.reason},'pending'))-1;
if ~isempty(pending)
    pool=gcp('nocreate');
    if isempty(pool)
        pool=parpool('Processes',min(n_workers,numel(pending)));
    end
    futures(1,numel(pending))=parallel.FevalFuture;
    for k=1:numel(pending)
        futures(k)=parfeval(pool,@local_run_pattern,1,pending(k),basep);
    end
    for k=1:numel(futures)
        [~,row]=fetchNext(futures);
        rows(row.pattern_id+1)=row;
        save(checkpoint,'rows','schema_version','signature','-v7.3');
        done=~strcmp({rows.reason},'pending');
        writetable(struct2table(rows(done)), ...
            fullfile(data_dir,'sign_redesign_patterns_partial.csv'));
        fprintf('sign %02d/63 %s ratio %.4f RMS %.4f peak %.1f runtime %.1fs\n', ...
            row.pattern_id,row.reason,row.max_ratio,row.rms_e,row.peak_u,row.runtime_s);
        drawnow;
        if ~row.accepted
            cancel(futures);
            error('Schema-30 sign scan stopped at failed pattern %d.',row.pattern_id);
        end
    end
end

T=struct2table(rows);
writetable(T,fullfile(data_dir,'sign_redesign_patterns.csv'));
save(fullfile(data_dir,'sign_redesign_patterns.mat'), ...
    'rows','schema_version','signature','-v7.3');
end

function row=local_run_pattern(code,basep)
signs=2*double(bitget(uint8(code),1:6))'-1;
p=basep;
p.sign_pattern=signs;
p.audit_detail=false;
p.save_stride=0.01;
p.max_wall_time=240;
timer=tic;
[o,~]=usv_simulate(p);
row=struct('pattern_id',code,'signs',sprintf('%+d',signs), ...
    'completed',o.success,'reason',o.reason,'T_final',o.T_final, ...
    'rms_e',o.rms_e_steady,'max_ratio',o.max_ratio,'peak_u',o.max_abs_u, ...
    'max_s',o.max_s,'max_theta',max(o.max_theta1,o.max_theta2), ...
    'envelope_event',o.envelope_event,'cap_triggered',o.cap_triggered, ...
    'state_guard',o.state_guard_triggered, ...
    'command_guard',o.command_guard_triggered, ...
    'ratio_guard',o.ratio_guard_triggered, ...
    'tracking_success',o.tracking_success, ...
    'margin_success',o.success && o.max_ratio<basep.max_ratio_guard, ...
    'accepted',o.tracking_success && o.success && o.max_ratio<basep.max_ratio_guard, ...
    'theorem_admissible',o.theorem_admissible,'runtime_s',toc(timer));
end
