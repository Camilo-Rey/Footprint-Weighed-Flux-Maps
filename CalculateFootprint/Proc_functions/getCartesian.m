function [LATsave,LONsave]=getCartesian(lat_gc,lon_gc,FX,FY)

Xsite = llh2ECEF(lat_gc,lon_gc,0);

% define the tranformation from the south east and up (SEZ) coordinate
% system to ECEF
t2 = -(pi/2+lon_gc);
t1 = -(pi/2-lat_gc);

R1 = [1,  0,       0;
    0,  cos(t1), sin(t1);
    0, -sin(t1), cos(t1)];

R3 = [ cos(t2), sin(t2), 0;
    -sin(t2), cos(t2), 0;
    0,       0,       1];

% Convert to lat long after merging
% rotation matrix
[FX2, FY2] = meshgrid(FX,FY);%  
[mtmp,ntmp] = size(FX2);
mat = R3*R1;

for j = 1:mtmp
        for k = 1:ntmp
            dist = [FX2(j,k);FY2(j,k)];            
            az1 = atan2(dist(2), dist(1));
            rho = norm(dist);            
            %dist_save = [dist_save,dist];
            
            % transform to the ECEF coordinate system
            rho_SEZ   = rho*[cos(az1);sin(az1);0];
            rho_ECEF  = mat*rho_SEZ;            
            Xtmp = Xsite + rho_ECEF;
            [lat,lon,h_ellp] = ECEF2llh(Xtmp);
            
            LATsave(j,k) = lat*180/pi;
            LONsave(j,k) = lon*180/pi;
        end              
end
end