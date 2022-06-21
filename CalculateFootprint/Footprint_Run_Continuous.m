%% Main Footprint_Run
% By Camilo Rey. Jul 2020.


%% Adjustments

if opt.DN==0;DT='night';
elseif opt.DN==1; DT='day';end  
if opt.BOTH==1; DT='both';end

Fall='';
for i=1:length(opt.FPcontour)
Flab=num2str(opt.FPcontour(i));
Fall=strcat(Fall,'-',Flab);
end

lega=Fall;% String of the customized contours to label files    
%% Load data

addpath('Proc_functions/')
addpath('Create_Site_Input/')

    
load(opt.L4file,'data','BADM')

data.WD=data.WD-opt.WDoffset; data.WD(data.WD<0)=data.WD(data.WD<0)+360;


% Load map
if opt.HaveMap==1
map_path=['..\Footprint_Map\' opt.site '\LandCover_' opt.site];
end
    
%% Index desired time

% Daytime vector
dayvec=data.ze<90;% daytime= 1. Nighttime=0

start_date=datenum(opt.start,'yyyy-mm-dd HH:MM');
end_date=datenum(opt.end2,'yyyy-mm-dd HH:MM');

% wind direction index (adjust if necessary)
wix=data.WD>opt.WDlow & data.WD<opt.WDhigh;

% Soil temp index if desired
timeT=data.time>opt.Timelow & data.time<opt.Timehigh;  

% Soil temp index if desired
%temp=data.ST>opt.STlow & data.ST<opt.SThigh;  

% % Stability index if desired
% data.zL=data.z./data.L;
% zL=data.zL>opt.zLlow & data.zL<opt.zLhigh;  


% FIX TEMP !!!!!!!!!!!!!!!!!!!!!!!

% Index all time within this period
if opt.BOTH==1
    ix=data.Mdate>=start_date & data.Mdate<=end_date & wix & timeT;
elseif opt.BOTH==0
    ix=data.Mdate>=start_date & data.Mdate<=end_date & dayvec==opt.DN & wix & timeT;
    ix2=data.Mdate>=start_date & data.Mdate<=end_date & dayvec==1 & wix & timeT;%daytime index
end

if exist('ix2','var') == 0; ix2=ix;end


   
%% Apply index to input variables:

UST=0.2;

    Mdate=data.Mdate(ix);
    plottime=datetime(datevec(Mdate));
    daynight=dayvec(ix);
    startH=data.time(ix)-BADM.window;
    endH=data.time(ix);
    try h=data.PBL_MS_gf(ix); catch h=ones(size(daynight))*750;end
    ustar=data.ustar(ix);
    ubar=data.ubar(ix);
    mbar=data.mbar(ix);
    vv=data.vv(ix);
    Lo=data.L(ix);Lo(ustar<UST)=nan;
    windir=data.WD(ix);
    DOY=data.DOY(ix);  
    year=data.year(ix);
    TA=data.TA(ix);
    ST=data.ST(ix);

    %Optional
    LE=data.LE(ix);LE(ustar<UST)=nan;
    FCH4=data.wm(ix);FCH4(ustar<UST)=nan;
    FCO2=data.wc(ix);FCO2(ustar<UST)=nan;
    %GPP=data.GPP(ix);GPP(ustar<UST)=nan;
    %if strcmp(site,'TZ'); GPP=data.GPP(ix);else;GPP=data.gpp_ANNnight(ix);end;GPP(ustar<UST)=nan;
    FH2O=data.wq(ix);FH2O(ustar<UST)=nan;
    Fuv=(data.ustar(ix)).^2;Fuv(ustar<UST)=nan;
    SH=data.H(ix);SH(ustar<UST)=nan;
        
    try
    GEB=(data.H(ix)+data.LE(ix))-data.RNET(ix); % Energy Balance
    LUE=GPP./data.PPFD(ix);% Light use efficiency
    WUE=GPP./data.wq(ix); % Water use efficiency
    catch
    GEB=nan(size(LE));
    LUE=nan(size(LE));
    WUE=nan(size(LE));
    end
    
    if strcmp(opt.site,'SWiso') 
    
    GEB=data.FC13(ix); 
    LUE=data.FC13_dsk(ix);
    WUE=data.FCO2iso(ix);  
    Fuv=data.FCO2iso_dsk(ix);
    ST=data.FC13_dsk(ix);
    end
    
    if strcmp(opt.site,'SWiso2') 
    
    GEB=data.ISO13C_RAC(ix); 
    LUE=data.ISO13C_RAC_dsk(ix);
    WUE=data.ISO13C_RAC_dsk2(ix);  
    Fuv=data.ISO18O_RAC_dsk(ix);
    ST=data.ISO18O_RAC_dsk(ix);
    end
    
        
    Canopy_height=data.z_veg(ix);
    Tower_height=data.z(ix);
    
    Mdate_day=data.Mdate(ix2);
    ustar_day=data.ustar(ix2);
    Lo_day=data.L(ix2);
    ubar_day=data.ubar(ix2);
    Tower_height_day=data.z(ix2);
    
    

%% load map or Define Grid

zm=nanmean(Tower_height)-nanmean(BADM.Canopy_height); % Average tower height
F2H=round(opt.f2h*zm,-2);

if opt.HaveMap==1
    load(map_path)
    PatchMap=MapB;
    
elseif opt.mscale==0
    
    lgt=F2H*2 ;% Lenght of one side of the square with the tower in the center
    PxSize=lgt/100;% Recommended 100*100 pixels
    FX=-lgt/2:PxSize:lgt/2;
    FY=-lgt/2:PxSize:lgt/2;
    PatchMap=nan;
    PatchCode='NA';
    Cmap=[1,1,1];
    
elseif opt.mscale==1
    
    lgt=F2H*2 ;% Lenght of one side of the square with the tower in the center    
    FX=[flip(-exp(0:0.1:log(lgt/2))),0, exp(0:0.1:log(lgt/2))];
    FY=[flip(-exp(0:0.1:log(lgt/2))),0, exp(0:0.1:log(lgt/2))];
    PxSize=round(nanmean(FX(2:end)-FX(1:end-1)),0);% Recommended no more than 5 m (m)
    PatchMap=nan;
    PatchCode='NA';
    Cmap=[1,1,1];
end

 %% Roughness

z=Tower_height;%   
z2=Tower_height_day;%   
d = 0.66*Canopy_height;
zo=0.1*Canopy_height;


%% Plotting adjustments
if opt.plotYN==1
PlotVec=ones(length(dayvec(ix)),1);% Ones for plot, zeros for no plot
else
PlotVec=zeros(length(dayvec(ix)),1);% Ones for plot, zeros for no plot
end            
    % Set a point of interest
    Tdist=ones(size(windir))*opt.distP; % 
    Twindir=ones(size(windir))*opt.dirP; % 
    

%% Run Footprint Code

    [COUNT,COMP,PATCHSUM_H,PATCHSUM_K,PATCHSUM_KM]=...
    FP_process_Aggregate(ustar,vv,Lo,windir,zo,z,d,ubar,FX,FY,PatchMap,Cmap,PlotVec,opt.site,DOY,year,h,startH,endH,...
    opt.FPcontour,opt.models,opt.mscale,Tdist,Twindir,FCH4,FCO2,FH2O,Fuv,SH,ST,GEB,LUE,WUE);

%% Variable Description:

metadata.variables = {'COUNT','Number of footprints in current subset','';...
             'PATCHSUM','A matrix with the percent contributions of each patch to the footprint (see PatchCode for correspondance)','';...
             'Mdate','Matlab date & time @ end of flux averaging period','';...
             'PatchCode','Code of each land cover within the opt.site (if available)','';...
             'footCount',' A matrix that counts the times that each cell was included in the footprint','';...
             'perF','Percent Footprint calculated from the cumulative FP','';...
             'MapAvg','Equals percF*100','';...
             'MapB','The original Map of patch distribution','';...
             'Cmap','Colormap used in the Map of the opt.site','';...
             'FCH4','CH4 flux','nmonl m-2 s-1';...
             'PxSize','The size of the pixel in the map used','m'};
         
%% Plot Climatology for current subset and save
addpath(genpath('googleearth'))    
% Compute coordinates to plot    

    global R_e Earth_E2 omega_e

    % radius of the Earth
    R_e = 6.37813649e6; % [m]

    % Earth's shape - eccentricity^2
    Earth_E2 = 0.006694385000^2;      % [-]

    % Mean Angular Velocity of the Earth
    omega_e =  7.29211585530e-5;    % [rad/s]
    
    lat_gc = BADM.lat_site*pi/180;
    lon_gc = BADM.lon_site*pi/180;
    [LATsave,LONsave]=getCartesian(lat_gc,lon_gc,FX,FY);
    [FX2, FY2] = meshgrid(FX,FY);%  
    RT= length(FY);
    
    
modLab={'Hsieh','Kljun','KM'};
    
for i=1:3
    FileLabel=[opt.site '-' modLab{i} '-' opt.start(1:4) opt.start(9:11) '-' opt.start(13:14) opt.start(16:17) '-to-' opt.end2(1:4) opt.end2(9:11) '-' opt.end2(13:14) opt.end2(16:17) '-' DT lega];
    if opt.models(i)==1
     
    % Plot Base Map
    figure;
    if ~isnan(PatchMap)
        pcolor(FX2,FY2(RT:-1:1,:),PatchMap*100);shading('flat');% This is the same as using flidud, which flips the rows so that the first one is the last one
        colormap(Cmap);
    else
        pcolor(FX2,FY2(RT:-1:1,:),nan(length(FX),length(FY)));shading('flat');% This is the same as using flidud, which flips the rows so that the first one is the last one
    end    
        hold on; xlabel('x [m]');  ylabel('y [m]')
        
    % Cumulative Footprint
    FG(i).percF=CalcPercF_fast(COMP(i).footCUM/COUNT)*100;
    [~,cA]=contour(FX2,FY2(RT:-1:1,:),FG(i).percF,opt.FPcontour,'Fill','off','ShowText','off');%shading('interp') 
    cA.LineColor = 'k'; cA.LineWidth = 1;hold on;   
    title(['Climatology ' FileLabel '%'])
    end
end
    legend(modLab)
for i=1:3
    FileLabel=[opt.site '-' modLab{i} '-' opt.start(1:4) opt.start(9:11) '-' opt.start(13:14) opt.start(16:17) '-to-' opt.end2(1:4) opt.end2(9:11) '-' opt.end2(13:14) opt.end2(16:17) '-' DT lega];
    % Export to kml
    if opt.models(i)==1
    figure
    kml_contour(LONsave,LATsave(RT:-1:1,:),FG(i).percF,...
    [opt.SaveDir 'kml\' FileLabel 'percent.kml'],opt.FPcontour,'r')
    end
end

%% Save all output    

clear data

    save([opt.SaveDir FileLabel opt.Sufix '.mat'])



