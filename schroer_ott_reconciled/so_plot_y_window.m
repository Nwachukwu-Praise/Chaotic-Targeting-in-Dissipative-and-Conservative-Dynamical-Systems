function yWindow = so_plot_y_window(cfg)
%SO_PLOT_Y_WINDOW Common y-range for every phase-space figure of a case.
pad = cfg.background.yPadding;
yMin = min(cfg.sourceRectangle.yMin, cfg.targetRectangle.yMin) - pad;
yMax = max(cfg.sourceRectangle.yMax, cfg.targetRectangle.yMax) + pad;
yWindow = [yMin, yMax];
end
