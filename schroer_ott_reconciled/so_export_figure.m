function paths = so_export_figure(fig, cfg, baseName)
%SO_EXPORT_FIGURE Write a figure in every requested format.
%
% cfg.figureFormats defaults to {'png'}.  Vector formats are useful for a
% thesis but the phase-portrait background is hundreds of thousands of
% points, so an SVG or EPS of these figures runs to tens of megabytes --
% add 'svg' or 'eps' deliberately, not by default.
formats = {'png'};
if isfield(cfg, 'figureFormats') && ~isempty(cfg.figureFormats)
    formats = cfg.figureFormats;
end
if ~exist(cfg.figureDirectory, 'dir')
    mkdir(cfg.figureDirectory);
end
paths = strings(0, 1);
for i = 1:numel(formats)
    fmt = lower(char(formats{i}));
    target = fullfile(cfg.figureDirectory, [char(baseName), '.', fmt]);
    switch fmt
        case {'png', 'jpg', 'tif', 'pdf'}
            exportgraphics(fig, target, 'Resolution', 200);
        case 'svg'
            exportgraphics(fig, target, 'ContentType', 'vector');
        case 'eps'
            print(fig, target, '-depsc', '-painters');
        otherwise
            warning('SchroerOtt:UnknownFigureFormat', ...
                'Skipping unknown figure format "%s".', fmt);
            continue;
    end
    paths(end + 1, 1) = string(target); %#ok<AGROW>
end
end
