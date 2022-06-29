function X_ECEF = llh2ECEF(lat,long,h_ellp)
%% DESCRIPTION
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%       Last updated:         April 21, 2011
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Given a position in geodetic latitude, longitude, and height above an
% ellipsoid determine the position in the Earth Centered Earth Fixed
% coordinate frame.
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
%    
%            lat = geodetic latitude                             [rad]
%           long = longitude                                     [rad]
%         h_ellp = height above an ellipsoidal model of Earth   *[length]
%
% -------------------------------------------------------------------------
% OUPUT
% -------------------------------------------------------------------------
%
%       X = Position vector of the spacecraft expressed in      *[length]
%           the ECEF coordinate frame.  
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   
% * this quantity can be expressed in either m or km or etc as long
%   as the global value of R_e (the Earth's Radius) is in consitant units.
%
%% DEFINE GLOBAL VARIABLES TO BE USED

global R_e Earth_E2

%% IMPLEMENTATION

N_phi  = R_e / sqrt(1 - Earth_E2*sin(lat)^2);

X_ECEF(1,1) = (N_phi+h_ellp)*cos(lat)*cos(long);
X_ECEF(2,1) = (N_phi+h_ellp)*cos(lat)*sin(long);
X_ECEF(3,1) = (N_phi*(1-Earth_E2)+h_ellp)*sin(lat);

