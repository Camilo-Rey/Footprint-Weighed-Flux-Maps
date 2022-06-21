%% Input footprint options

% By Camilo Rey. Jul 2020.

%%%%%% Start loading Footprint options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except data
       
    % What type of data?
    opt.DataType=1;% 1: biometeLab Mat files (ready), 2: Ameriflux CSV (pending)
    opt.site='SWiso2';
    opt.L4file=['..\L4_Processing\Processed\' opt.site '_L4.mat'];

    
    % how?
    opt.models=[0 0 1]; % Select any combination of 3 models [ Hsieh  Kljun  K&M ];
    opt.mscale=0; % The scale of the initial unrotated matrix. Must be zero for now
    opt.mscale2=0;% The scale of the final, rotated matrix. Must be zero for now

    opt.PBLdata=0; % 0 if there is no PBL data available
    opt.calcZo=1; %Calculate roughness length
    opt.calcH=1; % Calculate aerodynamic canopy heigth

    opt.HaveMap=0;% 1= Yes, I have a raster with land covers, 0= I dont have one
    opt.plotYN=0;% Ones for plot, zeros for no plot. % Normally stays at zero if running more than ~30 footprints

    % when?
    opt.BOTH=0; % Both daytime and nighttime at the same time, 1=yes, 0= no,separetely
    opt.DN=0;% daytime= 1. Nighttime=0; (only works if opt.BOTH==0)
    opt.PlotVec=0; % Plot each half hourly?
    
    % Defining an index to run the footprint in day of year        
    opt.start='2018-01-138 00:00'; % format YYYY-01-DOY HH:MM
    opt.end2='2018-01-154 24:00';
    
 
    opt.FPcontour=[50,80];% Choose Footprint Percentages. no more than 60-70 % for nighttime
    opt.f2h=100; % Fetch to height ratio defines the initial size of grid. varies from 100-150
    opt.pxS=200; % Pixels per side
    % Optional: Wind direction adjustment (if needed)
    opt.WDoffset=0;
    
    % Optional: Select a point around the tower to evaluate
    opt.distP=200; % Distance from the tower
    opt.dirP=360; % Azimuth
    
    % Optional: Filter based on wind direction or temperature
    opt.WDlow=0;
    opt.WDhigh=360;
    opt.STlow=0;
    opt.SThigh=35;   
    opt.Timelow=0;
    opt.Timehigh=24;
    
    % Saving options
    opt.SaveDir=(['..\Footprint_Output\' opt.site '\continuous\']);
    opt.Sufix='';
        
disp('Options recorded. Run Footprint_Run_half_hourly.m')  

%%%%%%%%% Stop entering footprint options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 