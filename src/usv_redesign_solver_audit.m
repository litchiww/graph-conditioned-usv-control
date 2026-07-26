function T=usv_redesign_solver_audit(root_dir)

if nargin<1, root_dir=fileparts(fileparts(mfilename('fullpath'))); end
cfg=usv_config(); data_dir=fullfile(root_dir,'data');
if ~exist(data_dir,'dir'), mkdir(data_dir); end

codes=[2 36 57 63];
strides=[0.01 0.05];
max_steps=[0.005 0.0025];
specs=[];
for code=codes
  for max_step=max_steps
    for stride=strides
        specs=[specs;code,stride,max_step];
    end
  end
end

basep=usv_redesign_candidate(cfg.get_params('C2_full'));
out_file=fullfile(data_dir,'sign_redesign_solver_audit.csv');
rows=table();
if exist(out_file,'file')
    previous=readtable(out_file,'TextType','string');
    keep=previous.completed & previous.reason=="none" & previous.T_final>=basep.T-1e-9;
    rows=previous(keep,:);
end
pending=true(size(specs,1),1);
for k=1:size(specs,1)
    pending(k)=~any(rows.pattern_id==specs(k,1) & ...
        abs(rows.save_stride-specs(k,2))<1e-12 & ...
        abs(rows.max_step-specs(k,3))<1e-12);
end
specs=specs(pending,:);
if ~isempty(specs) && license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if isempty(pool), pool=parpool('Processes',min(4,size(specs,1))); end
    futures(1,size(specs,1))=parallel.FevalFuture;
    for k=1:size(specs,1)
        futures(k)=parfeval(pool,@local_run_solver,1,basep,specs(k,:));
    end
    for k=1:numel(futures)
        [~,row]=fetchNext(futures);
        rows=[rows;row];
        writetable(rows,out_file);
        local_report(row);
    end
elseif ~isempty(specs)
    for k=1:size(specs,1)
        row=local_run_solver(basep,specs(k,:));
        rows=[rows;row];
        writetable(rows,out_file);
        local_report(row);
    end
end
T=sortrows(rows,{'pattern_id','max_step','save_stride'});
writetable(T,out_file);
end

function row=local_run_solver(basep,spec)
code=spec(1); stride=spec(2); max_step=spec(3);
p=basep;
p.sign_pattern=2*double(bitget(uint8(code),1:6))'-1;
p.audit_detail=false; p.save_stride=stride; p.max_step=max_step;
p.max_wall_time=900;
timer=tic; [o,~]=usv_simulate(p); runtime=toc(timer);
row=table(code,stride,max_step,o.success,string(o.reason),o.T_final, ...
    o.rms_e_steady,o.max_ratio,o.max_abs_u,o.max_s,runtime, ...
    'VariableNames',{'pattern_id','save_stride','max_step','completed','reason', ...
    'T_final','rms_e','max_ratio','peak_u','max_s','runtime_s'});
end

function local_report(row)
fprintf('sign %02d stride %.3g maxstep %.4g: %s ratio %.5f RMS %.5f peak %.1f\n', ...
    row.pattern_id,row.save_stride,row.max_step,row.reason, ...
    row.max_ratio,row.rms_e,row.peak_u);
drawnow;
end
