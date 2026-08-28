function profile = so_new_performance_profile()
%SO_NEW_PERFORMANCE_PROFILE Counters required by the first-light run.
profile.forwardFamilyBuilds = 0;
profile.backwardFamilyBuilds = 0;
profile.backwardFamilyCacheReuses = 0;
profile.indexedCrossingQueries = 0;
profile.exactSegmentTestsAfterSpatialFiltering = 0;
profile.naiveComparisonsAvoided = 0;
profile.jProbesPruned = 0;
profile.runtimeByStage = [];
profile.peakCurvePointCount = 0;
end

