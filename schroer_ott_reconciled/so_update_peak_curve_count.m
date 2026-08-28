function profile = so_update_peak_curve_count(profile, family)
%SO_UPDATE_PEAK_CURVE_COUNT Track largest adaptive curve encountered.
for i = 1:numel(family)
    if isfield(family{i}, 'pointCount')
        profile.peakCurvePointCount = max(profile.peakCurvePointCount, family{i}.pointCount);
    end
end
end

