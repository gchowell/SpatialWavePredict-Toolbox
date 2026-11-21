function [value,isterminal,direction]=spatialWaveActivationEvent(t,x,onset_thr,lastActive)
%SPATIALWAVEACTIVATIONEVENT Event: the last active sub-epidemic reaches the threshold.
%   Used by simulateSpatialWave.m. The integration stops when the cumulative
%   count of sub-epidemic LASTACTIVE rises through ONSET_THR; the simulator
%   then activates the next one and starts a new segment. Nothing is changed
%   in this function.
value=x(lastActive)-onset_thr;
isterminal=1;
direction=1; % upward crossings only
end
