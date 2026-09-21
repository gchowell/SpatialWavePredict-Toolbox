function datenum1=getSpatialWaveVerificationDate(caddate1,DT)
%GETSPATIALWAVEVERIFICATIONDATE Date argument used to load verification data.
%
% caddate1 is the date of the final observation used for calibration.
% For daily and weekly series, getData selects the observation represented
% by its date argument, so advance by one sampling interval. For annual
% series, getData selects the year after its date argument internally, so
% retain the calibration-end year to avoid skipping a year.

narginchk(2,2);

if isstring(caddate1)
    if ~isscalar(caddate1)
        error('SpatialWave:InvalidCalibrationDate', ...
            'caddate1 must be a scalar date in mm-dd-yyyy format.');
    end
    caddate1=char(caddate1);
end

if ~ischar(caddate1) || size(caddate1,1)~=1 || ...
        isempty(regexp(caddate1,'^\d{2}-\d{2}-\d{4}$','once'))
    error('SpatialWave:InvalidCalibrationDate', ...
        'caddate1 must be a scalar date in mm-dd-yyyy format.');
end

validateattributes(DT,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'DT',2);

try
    datenum1=datenum(caddate1,'mm-dd-yyyy');
catch
    error('SpatialWave:InvalidCalibrationDate', ...
        'caddate1 must be a valid date in mm-dd-yyyy format.');
end

if DT~=365
    datenum1=datenum1+DT;
end

end
