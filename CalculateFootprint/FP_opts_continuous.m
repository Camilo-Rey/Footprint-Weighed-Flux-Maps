%% Input footprint options

% By Camilo Rey-Sanchez. Last modified Jun 2022.
% Before running this file, an L4 file must be created using the code in
% the folder L4_processing.

%%%%%% Start loading Footprint options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except data
       
    % What type of data?
    opt.site='NC4';% Short name or acronym for the site of interest
    opt.L4file='E:\FluxData\L4_Processing\Processed\NC4_L4.mat';% The directory of the L4 file
    
    % how?
    opt.models=[0 0 1]; % Select any combination of 3 models [ Hsieh  Kljun  K&M ];
    opt.mscale=0; % The scale of the initial unrotated matrix. Must be zero for now
    opt.mscale2=0;% The scale of the final, rotated matrix. Must be zero for now

    opt.PBLdata=0; % 0 if there is no PBL data available
    opt.HaveMap=0;% Default: 0= I dont have one. 1= Yes, I have a raster with land covers, 
    opt.plotYN=0;% Ones for plot, zeros for no plot. % Normally stays at zero if running more than ~30 footprints

    % when?
    opt.BOTH=0; % Both daytime and nighttime at the same time, 1=yes, 0= no,separetely
    opt.DN=0;% daytime= 1. Nighttime=0; (only works if opt.BOTH==0)
    opt.PlotVec=0; % Plot each half hourly?
    
    % Defining an index to run the footprint in day of year        
    opt.start='2018-01-138 00:00'; % format YYYY-01-DOY HH:MM
    opt.end2='2018-01-144 24:00'; % format YYYY-01-DOY HH:MM    
 
    opt.FPcontour=[50,80];% Choose Footprint Percentages. no more than 60-70 % for nighttime
    opt.f2h=100; % Fetch to height ratio defines the initial size of grid. varies from 100-150
    opt.pxS=200; % Pixels per side of the square. default=200
    
    % Optional: Wind direction adjustment (if needed)
    opt.WDoffset=0;
    
    % Optional: Select a point around the tower to evaluate
    opt.distP=00; % Distance from the tower (m)
    opt.dirP=00; % Azimuth
    
    % Optional: Filter the dataset for footprint aggregation based on wind direction or temperature
    % Leave default values if no filter is to be applied.
       
    opt.WDlow=0;% Wind direction, default 0
    opt.WDhigh=360; % Wind direction high, default 360
    opt.Timelow=0;% Time in hours, default 0
    opt.Timehigh=24;% Time in hours, default 24
    
    %% Add fluxes for footprint-weighed flux maps
    disp('loading data. Please wait a few seconds...')
    load(opt.L4file)
    Fuv=(data.ustar).^2;% Momentum flux ( m-2 s-2)
    % Add fluxes here. If the default flux does not exist,
    % it can be left as is or replaced by other flux here:
    try FS0=data.wm;opt.FS0='CH4 flux (nmol m-2 s-1)';catch FS0=nan(size(Fuv));end % CH4 flux (nmol m-2 s-1)
    try FS1=data.wc;opt.FS1='CO2 flux (umol m-2 s-1)';catch FS1=nan(size(Fuv));end % CO2 flux (umol m-2 s-1)
    try FS2=data.wq;opt.FS2='H2O flux (mmol m-2 s-1)';catch FS2=nan(size(Fuv));end % H2O flux (mmol m-2 s-1)
    FS3=Fuv;opt.FS3='Momentum flux ( m2 s-2)'; % Momentum flux ( m2 s-2)
    try FS4=data.H;opt.FS4='Sensible heat flux (W m-2)'; catch FS4=nan(size(Fuv));end % Sensible heat flux (W m-2)
    try FS5=data.LE;opt.FS5='Latent heat flux (W m-2)';catch FS5=nan(size(Fuv));end % Latent heat flux (W m-2)
    FS6=nan(size(Fuv));opt.FS6='None'; % Add other flux here
    FS7=nan(size(Fuv));opt.FS7='None'; % Add other flux here
    FS8=nan(size(Fuv));opt.FS8='None'; % Add other flux here
    
    %% Saving options
    opt.SaveDir=(['E:\Footprint_Output\' opt.site '\continuous\']);
    opt.Sufix='test';% A sufix to add to the name of the file
        
disp('Options recorded. Run Footprint_Run_Continuous.m')  

%%%%%%%%% Stop entering footprint options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 