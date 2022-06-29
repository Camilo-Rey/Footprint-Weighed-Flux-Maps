    function [ZE]=CalculateZenith(lon_site,lat_site,altitude,Mdate,UTC_offset)
    
    Dtime=Mdate;

    location.longitude = lon_site; % degrees - longitude of tower location
    location.latitude  =  lat_site; % degrees - latitude of tower location
    location.altitude  =  altitude; % m - altitude of tower location    
    
        %compute local sun and zenith angle (refer time to Greenwich)
        gtG = -(UTC_offset/24); % Time offset to get to Greenwhich
        
        ZE=nan(length(Dtime),1);
        for i=1:length(Dtime)
        zx = datestr(Dtime(i)+gtG,'dd-mmm-yyyy HH:MM:SS');
        sun = sun_position(zx, location);
        ZE(i)=sun.zenith;
        end
        
    end