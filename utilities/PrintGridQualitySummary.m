function PrintGridQualitySummary(AreaRatioQuality, ShapeQuality, AngleQuality, Skewness)
disp('===============  网格质量统计  ===============');
%%
EquiAngleSkewness       = 1-Skewness(1);
averEquiAngleSkewness   = 1-Skewness(2);
disp(['maxEquiAngleSkewness      = ', num2str(EquiAngleSkewness, '%.3f')]);
disp(['averEquiAngleSkewness     = ', num2str(averEquiAngleSkewness, '%.3f')]);

%%
minQ  = ShapeQuality(1);
maxQ  = ShapeQuality(2);
averQ = ShapeQuality(3);
disp('---------------------------------------------');
%             disp(['maxQuality    = ', num2str(maxQ, '%.3f')]);
disp(['minQuality    = ', num2str(minQ, '%.3f')]);
disp(['averQuality   = ', num2str(averQ, '%.3f')]);

%%
minAngle = AngleQuality(1);
maxAngle = AngleQuality(2);
averMinAngle = AngleQuality(3);
averMaxAngle = AngleQuality(4);
disp('---------------------------------------------');
disp(['minAngle     = ', num2str(minAngle, '%.3f')]);
disp(['averMinAngle = ', num2str(averMinAngle, '%.3f')]);
disp('---------------------------------------------');
disp(['maxAngle     = ', num2str(maxAngle, '%.3f')]);
disp(['averMaxAngle = ', num2str(averMaxAngle, '%.3f')]);

%%
minAR  = AreaRatioQuality(1);
maxAR  = AreaRatioQuality(2);
averAR = AreaRatioQuality(3);
disp('---------------------------------------------');
%             disp(['maxAreaRatio  = ', num2str(maxAR, '%.3f')]);
disp(['minAreaRatio  = ', num2str(minAR, '%.3f')]);
disp(['averAreaRatio = ', num2str(averAR, '%.3f')]);
disp('---------------------------------------------');
disp(' ');
end