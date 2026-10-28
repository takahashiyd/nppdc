%% run_tutorial.m -- runs nppdc_tutorial.m and saves its three figures.
try
    root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
    addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);
    nppdc_tutorial
    figs = findobj('Type', 'figure');
    % order by creation
    [~, ord] = sort([figs.Number]);
    figs = figs(ord);
    names = {'ipdc_validation.png', 'ipdc_asymp.png', 'asymp_calibration.png'};
    for k = 1:min(numel(figs), numel(names))
        out = fullfile('figures', names{k});
        style_paper_figure(figs(k));
        set(figs(k), 'PaperPositionMode', 'auto');
        print(figs(k), out, '-dpng', '-r150');
        fprintf('SAVED_FIGURE: %s\n', out);
    end
catch ME
    fprintf(2, 'ERROR: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    exit(1);
end
