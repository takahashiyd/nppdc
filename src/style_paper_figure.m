function style_paper_figure(fig)
%STYLE_PAPER_FIGURE  Apply consistent publication-quality formatting.
%
%   style_paper_figure(fig)  styles all axes and legends in fig with:
%     - Helvetica font, size 12 (axes) / 10 (legend)
%     - Ticks out, top/right axis lines hidden, axes line width 0.8
%     - Soft grid (alpha ~0.12), grid behind data
%     - Title weight normal (not bold)
%     - Legend with no box, transparent background
%     - White figure background
%
%   Call this AFTER all plotting commands and BEFORE printing. It does not
%   touch line colors set explicitly in plot calls (use the colour codes
%   returned by PAPER_PALETTE for that).

    if nargin < 1 || isempty(fig), fig = gcf; end
    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');

    ax_list = findall(fig, 'Type', 'axes');
    for h = ax_list(:).'
        try
            set(h, ...
                'FontName',        'Helvetica', ...
                'FontSize',        12, ...
                'LineWidth',       0.8, ...
                'TickDir',         'out', ...
                'TickLength',      [0.015 0.025], ...
                'Box',             'off', ...
                'XGrid',           'on', ...
                'YGrid',           'on', ...
                'GridColor',       [0 0 0], ...
                'GridAlpha',       0.12, ...
                'MinorGridAlpha',  0.06, ...
                'Layer',           'top');
            t = get(h, 'Title');
            if ~isempty(t)
                set(t, 'FontWeight', 'normal', 'FontSize', 12);
            end
        catch
        end
    end

    legs = findall(fig, 'Type', 'Legend');
    for h = legs(:).'
        try
            set(h, ...
                'FontName', 'Helvetica', ...
                'FontSize', 10, ...
                'Box',      'off', ...
                'Color',    'none', ...
                'EdgeColor','none');
        catch
        end
    end

    sgts = findall(fig, 'Type', 'subplottext');
    for h = sgts(:).'
        try
            set(h, 'FontName', 'Helvetica', 'FontWeight', 'normal', 'FontSize', 13);
        catch
        end
    end

    % Force all text annotations (including tiledlayout titles set via
    % title(tl, ...)) to a normal-weight Helvetica.
    all_text = findall(fig, 'Type', 'text');
    for h = all_text(:).'
        try
            set(h, 'FontName', 'Helvetica', 'FontWeight', 'normal');
        catch
        end
    end
end
