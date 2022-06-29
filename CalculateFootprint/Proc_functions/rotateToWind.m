function [xr,yr] = rotateToWind(FX,FY,windir,F, plotYN)
%       windir - wind direction, from north, degrees
%       
%        FX - meshgrid, x coordinates, of mask map
%
%        FY - meshgrid, y coordinates, of mask map
%
%        F - mask map on FX,FY coordinates

if nargin==3
    plotYN=0;
end

wd = (windir -90)/180*pi;
xr = FX*cos(wd) + FY*sin(wd);
yr = FY*cos(wd) - FX*sin(wd);

if plotYN == 1 && nargin > 3
     figure(2)
    pcolor(xr,yr,F);shading('interp')
    xlabel('distance upwind [m]')
    ylabel('distance crosswind [m]')
    title('Area of interest')
end


end