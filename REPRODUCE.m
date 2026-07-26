function REPRODUCE(stage)
if nargin<1, stage='figures'; end
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'src'));
switch lower(stage)
    case 'verify'
        VERIFY_PACKAGE;
    case 'figures'
        MAIN('plots');
    case 'canonical'
        MAIN('math');
        MAIN('cases');
        MAIN('lambda');
        MAIN('hydro');
        MAIN('closedloop');
        MAIN('plots');
    case 'closedloop'
        MAIN('closedloop');
    case 'signs'
        if license('test','Distrib_Computing_Toolbox')
            usv_redesign_sign_parallel(root,4);
        else
            usv_validation_suite('signs_redesign',root);
        end
        usv_direction_identification_audit(root);
    case 'random'
        usv_redesign_mc_audit(root,30);
    case 'solver'
        usv_redesign_solver_audit(root);
    case 'long'
        usv_redesign_long_audit(root);
    case 'audits'
        REPRODUCE('signs');
        REPRODUCE('random');
        REPRODUCE('solver');
        REPRODUCE('long');
    case 'all'
        REPRODUCE('canonical');
        REPRODUCE('audits');
        MAIN('plots');
    otherwise
        error('Unknown reproduction stage: %s',stage);
end
end
