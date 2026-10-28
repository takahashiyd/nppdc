function run_all(which)
%RUN_ALL  Regenerate every figure and table of the tutorial into figures/.
%
%   run_all            runs all demos (several hours: the MC scripts dominate)
%   run_all({'bivar_vs_multivar','hawkes_demo'})   runs a subset
%
%   Figure / table -> script map (figures are written to figures/):
%     ipdc_validation.png, ipdc_asymp.png, asymp_calibration.png  run_tutorial (wraps nppdc_tutorial)
%     bivar_vs_multivar.png                                       bivar_vs_multivar
%     ma_bandpass.png                                             ma_bandpass_example
%     bias_correction.png                                         validate_bias_correction
%     check_power.png, check_nongauss.png, Table tab:alpha_sweep  more_checks
%       (also check_alpha_sweep.png, a plot of that table; not in the paper)
%     hawkes_validation.png                                       hawkes_demo
%     hawkes_calibration.png                                      hawkes_mc
%     hodge_hierarchy.png                                         hodge_demo  (needs data/allen)
%     allen_lfp.png                                               allen_lfp_analysis  (needs data/allen)
%     Table tab:stress (console)                                  stress_test
%       (also stress_A_scale_mismatch.png, stress_B_peaky.png; not in the paper)
%     preconditioning table, Sec. 4.3 (console)                   precond_residual
%     Table tab:bias_small_m (console)                            test_bias_small_m
%     DOF calibration, Sec. 9.5 (console)                         calibrate_dof
%     Carter bias check, Sec. 9.6 (console)                       verify_carter_mean
%     Hodge calibration, Sec. 8.3 (console)                       hodge_jack_validate
%
%   Each demo is a script that sets its own paths and cd's to this folder,
%   so they can also be run one at a time from the MATLAB editor.
%   extract_allen_subset rebuilds data/allen/subset_movie.mat from the raw
%   NWB files and is not run here.

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'), fullfile(root, 'demos'));
    cd(root);
    if ~exist('figures', 'dir'), mkdir('figures'); end

    all_demos = {'bivar_vs_multivar', 'run_tutorial', 'ma_bandpass_example', ...
                 'validate_bias_correction', 'more_checks', ...
                 'hawkes_demo', 'hawkes_mc', 'hodge_demo', 'stress_test', ...
                 'precond_residual', 'allen_lfp_analysis', 'test_bias_small_m', ...
                 'calibrate_dof', 'verify_carter_mean', 'hodge_jack_validate'};
    if nargin < 1 || isempty(which), which = all_demos; end
    if ischar(which), which = {which}; end

    logf = fopen(fullfile(root, 'figures', 'run_all.log'), 'a');
    fprintf(logf, '\n==== run_all %s ====\n', datestr(now));
    for k = 1:numel(which)
        name = which{k};
        fprintf('\n=== %s ===\n', name); fprintf(logf, '=== %s ===\n', name);
        t0 = tic;
        try
            evalin('base', name);
            fprintf(logf, 'OK   (%.0f s)\n', toc(t0));
        catch ME
            fprintf(2, '  FAILED: %s\n', ME.message);
            fprintf(logf, 'FAIL %s\n', ME.message);
        end
        close all
    end
    fclose(logf);
    fprintf('\n=== DONE ===\n');
end
