%% L4 Processing
% By Camilo Rey Sanchez
% Jan 27 2021
%% This code uses input Ameriflux .csv. data and adds extra information including:
% Fix tower height with respect to changes in water level (if available)
% Calculate aerodynamic canopy height
% Correct fluxes based on enery balance approach
% Add PBL data to the main data file
% convert Ameriflux .csv. to Matlab format and create a L4 file

clear all

SaveDir='E:\FluxData\L4_Processing\Processed\';

% addpath(genpath('Q:\My Drive\FluxProcess\eddycovariancegeneral'));
% addpath(genpath('Q:\My Drive\UC-Berkeley\Footprint_share\CalculateFootprint\Proc_functions'))
% addpath(genpath('Q:\My Drive\DATA\Subroutines'));

%% Enter options for the site

    site='NC4';

    cluster='other';
    BADM.filename='AMF_US-NC4_BASE_HH_5-5.csv';
    BADM.filepath='Q:\My Drive\NC-State\Methane_Alligator\data\';
    BADM.window=0.5;% Time BADM.window of the data in hours (0.5 for half-hourly)    
    BADM.lat_site = 35.787717;
    BADM.lon_site = -75.903952;
    BADM.altitude=0; % sea level
    BADM.UTC_offset=-5;
    BADM.Tower_height=35;% m
    BADM.Canopy_height=17.5;% m  

    % Read data
    disp('Reading Ameriflux csv standard. Please wait a few seconds')
    [AF]=Read_AF_StandFormat(BADM.filename,BADM.filepath);      
    % Get zenith for daytime vector)  
    disp('Calculating Zenith angles. This may take a couple of minutes')
    [ze]=CalculateZenith(BADM.lon_site,BADM.lat_site,BADM.altitude,AF.Mdate,BADM.UTC_offset);        
    data.ze=ze;
    data.time=hour(AF.dateplot)+(minute(AF.dateplot)/60);           
    data.Mdate=AF.Mdate;    
    
    % Update Labels if needed (e.g. from 1_1_1 to 1_1_2)
    data.LE=AF.LE_1_1_1;
    data.H=AF.H_1_1_1;
    data.TA=AF.TA_1_1_1;
    data.ST=AF.TS_1_1_1;
    data.RH=AF.RH_1_1_1;
    data.VPD=AF.VPD_PI_1_1_1;
    data.G=AF.G_1_1_1;
    data.RNET=AF.NETRAD_1_1_1;   
    data.vv = AF.V_SIGMA_1_1_1.^2;
    data.ustar = AF.USTAR_1_1_1;   % 
    data.wc = AF.FC_1_1_1;   %       
    data.PPFD=AF.PPFD_IN_1_1_1;
    Lvap = (2.501-0.00237*data.TA)*1E6; % J/kg latent heat of vaporiz. - Eqn in back of Stull pg. 641   
    wq_mmol=(data.LE./Lvap)*1000*1000/18.0153;% mmol m-2 s-1
    param.Lt= 2260; % Latent heat of vaporization j/g
    wq_mmol2=(data.LE/param.Lt)*1000/18.0153;% mmol m-2 s-1
    data.wq = wq_mmol;        
    data.WL=AF.WTD_1_1_1; % in m
    data.wm = AF.FCH4_1_1_1;   
    data.mbar= nan(size(AF.Mdate)); 
    data.ubar=AF.WS_1_1_1;    
    data.L=BADM.Tower_height./AF.ZL_1_1_1;    
    data.WD=AF.WD_1_1_1;    
    data.zL=AF.ZL_1_1_1; 
    data.DOY=day(datetime(datevec(data.Mdate)),'dayofyear');
    data.year=year(datetime(datevec(data.Mdate)));    

    plottime=datetime(datevec(data.Mdate));
    % Daytime vector
    dayvec=data.ze<90;% daytime= 1. Nighttime=0
    data.daynight=data.ze<90;   


%% Try water level (update if available)

% data.z=BADM.Tower_height-data.WL;

data.z=ones(size(data.Mdate))*BADM.Tower_height;   

%% Aerodynamic canopy height from PeenyPacker and Baldocchi
daysAVG=30;%days to average
[h_dsk,ZvegSmooth]=CalculateAerodynamicCanopyHeight(data.Mdate,data.ustar,data.z,data.L,data.ubar,daysAVG,1);

%ZvegSmooth(~isnan(h_dsk))=h_dsk(~isnan(h_dsk));
ZvegSmooth(ZvegSmooth<0)=nan;
data.z_veg=ZvegSmooth;

data.zo=0.1*data.z_veg;% Roughness length
data.d0 = 0.66*data.z_veg; % m - zero plane dispacement height
data.z_abv = (BADM.Tower_height-data.z_veg); % m - tower height above veg height (indicated in opts structure). Note, this is the aerodynamic veg height, not actual, and is used to create non-dimensional cospectra. 

[d3,zo3,n3,stat3]=RoughnessLength(BADM.Tower_height,BADM.Canopy_height,data.L,data.ustar,data.ubar,data.daynight);
data.zo3=ones(length(data.ubar),1)*zo3;  % m - roughness length on Mauer et al (2013) from Gil Bohrer's lab.


%% Add PBL data
% (update if available)
        
    data.h=nan(size(data.Mdate));
    data.h_gf=nan(size(data.Mdate));    
    

%% Wind Roses to check conditions

    if sum(~isnan(data.ubar))>5
    WindRose(data.WD,data.ubar,'anglenorth',0,'angleeast',90,'labels',{'N (0°)','S (180°)','E (90°)','W (270°)'},'freqlabelangle',45,...
        'titlestring', {'Wind rose current subset (without offset)';' '});
    else
        disp('Too little data to plot wind rose')
    end


%% Save

save([SaveDir site '_L4.mat'],'-v7.3','data','BADM')

    