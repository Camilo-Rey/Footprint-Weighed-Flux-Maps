function [xout,nspike] = despike(x,n,lim)

% DESPIKE Basic despiking using a median filter.
% Spikes are detected if the difference of the data point to the filtered
% time series exceeds a chosen limit.
%
% INPUTS:
% x:   data vector with spikes
% n:   order of median filter (n=3 detects spikes consisting of
%      single data points only)
% lim: detection limit (i.e. spike if |x(i)-x_filtered(i)|>lim)
%
% OUTPUTS:
% xout: despiked data vector
% nspike: number of detected spikes
%
% Required MATLAB toolboxes: Signal Processing Toolbox

% AUTHOR: Patrick Sturm <pasturm@ethz.ch>
%
% COPYRIGHT 2011 Patrick Sturm
% This file is part of Eddycalc.
% Eddycalc is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% For a copy of the GNU General Public License, see
% <http://www.gnu.org/licenses/>.

% extend data vector by n data points at both ends
xi = [flipud(x(1:n)); x; flipud(x(end-n+1:end))];

% apply order n median filter
xfilt = medfilt1(xi,n);

% difference to filterd data
dif = abs(xfilt(1+n:end-n)-x);

% number of spikes
nspike = sum(dif>lim);

% replace spikes with NaNs
x(dif>lim) = NaN;
xout = x;