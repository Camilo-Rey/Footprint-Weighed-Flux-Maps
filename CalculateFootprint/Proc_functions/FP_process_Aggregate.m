function [ncalc,COMP,patchsum_H,patchsum_K,patchsum_KM]=...
    FP_process_Aggregate(ustar,vv,Lo,windir,zo,z,d,ubar,FX,FY,Map,Cmap,PlotVec,site,DOY,year,h,startH,endH,...
    contourMarks,models,mscale,Tdist,Twindir,FCH4,Fs1,Fs2,Fs3,Fs4,Fs5,Fs6,Fs7,Fs8)
%% Input:
    % ustar: friction velocity (m s-1)
    % vv: variance of the vertical wind component. Needed to calculate
    %   sigma v for lateral dispersion of the footprint along the x axis
    % Lo= Obukhov length
    % windir= wind direction data
    % zo= roughness length
    % d= displacment height
    % ubar= main rotated wind speed
    % FX, FY, Map= X and Y coordinates of the grid for which to calculate
    %   the footprint
    % Map: A raster with codes pertaining to different land covers in the
    %   grid
    % Cmap: The conventions for the codes used in Map
    % PlotVec: A vector of zeros and ones indicating which half hours in
    %   the current subset should be plotted
    % Site, DOY, Year
    % h: Planetaty Boundary Layer height (m)
    % u: wind speed 
    % StartH: Start date and time of the first footprint
    % contourMarks: Percentage footprints to draw. The last value is use to
    % crop the final footprint.
    % models: a vector of three options for which models to run:
    %   [1, 1, 1] --> [Hsieh, Kljun, K&M]
    % mscale: not used yet
    % Tdist: The distance of one point of interest to calculate footprint
    %   values for
    % Twindir: The wind direction of one point of interest to calculate footprint
    %   values for    
    % endH: End time of the last footprint
    % FCH4: Half-hourly methane flux (umol m-2 s-1)
    % FsX= Type of flux to be weighted (CO2, H2O, H, LE, UV...), where X is the flux number from 1
    %   the number of n species to be tested.


%% Output:

    % ncalc: Number of footprints that were calculated
    % COMP: Structure of compiled, aggregated footprint results
        % patchsum: Sum of the source function for a given patch. 
        %   cc
        % footCUM: Sum of the source function in each cell for all the half hour (_H
        %  is for Hsieh and _K is for Kjlun.     
        % footCount: the number of times each tile is counted(Number of tiles in the footprint that
        %   are not nan, usually the 80% FP)
        % comp: Average CH4 flux distribution. The sum of the lumped CH4 for
        %   each tile divided by the number of tiles each tile was counted


sv=sqrt(vv);

%% Pre-allocating

patchsum_H=nan(length(ustar),length(PatchIndex(Map)));
patchsum_K=nan(length(ustar),length(PatchIndex(Map)));
patchsum_KM=nan(length(ustar),length(PatchIndex(Map)));

COMP=[];

% Name of the model
COMP(1).name='H';
COMP(2).name='K';
COMP(3).name='KM';

for i=1:length(models)
COMP(i).xmax=nan(length(ustar),1);% Distance of the maximum footprint contribution
COMP(i).footCUM=nan(length(FX),length(FY));% Cumulative footprint values

COMP(i).Fs0tile=zeros(length(FX),length(FY));% Flux map footprint-weighed
COMP(i).Fs1tile=zeros(length(FX),length(FY));
COMP(i).Fs2tile=zeros(length(FX),length(FY));
COMP(i).Fs3tile=zeros(length(FX),length(FY));
COMP(i).Fs4tile=zeros(length(FX),length(FY));
COMP(i).Fs5tile=zeros(length(FX),length(FY));
COMP(i).Fs6tile=zeros(length(FX),length(FY));
COMP(i).Fs7tile=zeros(length(FX),length(FY));
COMP(i).Fs8tile=zeros(length(FX),length(FY));
COMP(i).Zvegtile=zeros(length(FX),length(FY));

COMP(i).Fs0tile_abs=nan(length(FX),length(FY));% Flux map for CH4 absolute values
COMP(i).Fs1tile_abs=nan(length(FX),length(FY));
COMP(i).Fs2tile_abs=nan(length(FX),length(FY));
COMP(i).Fs3tile_abs=nan(length(FX),length(FY));

COMP(i).footCount=zeros(length(FX),length(FY));% number of time a cell within the end countourMark is counted
COMP(i).footCountS0=zeros(length(FX),length(FY));% number of time a cell within the end countourMark is counted for CH4 footprints
COMP(i).patchsum=nan(length(ustar),length(PatchIndex(Map)));% Percent of each patch in the footprint
COMP(i).pxR=nan(length(ustar),1);% Footprint value for the cell of interest
COMP(i).footWeight=zeros(length(FX),length(FY));
COMP(i).footWeightS0=zeros(length(FX),length(FY));
end


long=length(ustar);
ncalc=0;


%% Loop for every half hour (or any other time stamp)

%Lo(abs(Lo)<0.2)=nan;% avoiding highly stable conditions

for i=1:long 
if ~isnan(ustar(i)*Lo(i)*sv(i)*zo(i)*z(i)*d(i)*windir(i)*h(i)*ubar(i))
    ncalc=ncalc+1;
    disp(['iteration # ' num2str(i), ', Footprint # ' num2str(ncalc)])   
    
    
    % Comparison between Hsieh and Kjlun
    [FPall,pxR]=...
       Calc_FP_basic(ustar(i),Lo(i),sv(i),zo(i),z(i),d(i),windir(i),PlotVec(i), FX, FY, Map,Cmap, h(i),...
       contourMarks,Tdist(i),Twindir(i),ubar(i),models,mscale);
    
    if PlotVec(i)==1
    title([num2str(site) ', Hour: ' num2str(startH(i)) ' to ' num2str(endH(i)) ', DOY ' num2str(DOY(i)) ', year ' num2str(year(i))  ', u*=' num2str(round(ustar(i),2)), ', L=' num2str(round(Lo(i),2))])
    end
    
    if models(1)==1;COMP(1).pxR(i)=pxR.H;end
    if models(2)==1;COMP(2).pxR(i)=pxR.K;end
    if models(3)==1;COMP(3).pxR(i)=pxR.KM;end
    
    
    for j=1:3
        if models(j)==1
        COMP(j).xmax(i)=FPall(j).Xmax;
        tmp = cat(3,COMP(j).footCUM,FPall(j).foot);
        COMP(j).footCUM = nansum(tmp,3);
        FPall(j).fCrop(FPall(j).fCrop==0)=nan;
        
        try
        COMP(j).patchsum(i,:)=FPall(j).PatchSum; 
        catch
        end
        
        %% Create footprin-weighed flux map
        
        CH4_gen=(FPall(j).fCrop)*FCH4(i);
        Fs1_gen=(FPall(j).fCrop)*Fs1(i);
        Fs2_gen=(FPall(j).fCrop)*Fs2(i);
        Fs3_gen=(FPall(j).fCrop)*Fs3(i);
        Fs4_gen=(FPall(j).fCrop)*Fs4(i);
        Fs5_gen=(FPall(j).fCrop)*Fs5(i);        
        Fs6_gen=(FPall(j).fCrop)*Fs6(i);
        Fs7_gen=(FPall(j).fCrop)*Fs7(i);
        Fs8_gen=(FPall(j).fCrop)*Fs8(i);
        Zveg_gen=(FPall(j).fCrop)*d(i)/0.66;
        
        % Lumped scalar flux. Put scalar flux in each tile where the 80% FP is located
        CH4_genL=~isnan(FPall(j).foot)*FCH4(i);
        s1_genL=~isnan(FPall(j).foot)*Fs1(i);
        s2_genL=~isnan(FPall(j).foot)*Fs2(i);
        s3_genL=~isnan(FPall(j).foot)*Fs3(i);

        %% aggregate the flux maps
        tmpA = cat(3,COMP(j).Fs0tile,CH4_gen); COMP(j).Fs0tile = nansum(tmpA,3);
        tmpB = cat(3,COMP(j).Fs1tile,Fs1_gen); COMP(j).Fs1tile = nansum(tmpB,3);
        tmpC = cat(3,COMP(j).Fs2tile,Fs2_gen); COMP(j).Fs2tile = nansum(tmpC,3);
        tmpD = cat(3,COMP(j).Fs3tile,Fs3_gen); COMP(j).Fs3tile = nansum(tmpD,3);
        tmpE = cat(3,COMP(j).Fs4tile,Fs4_gen); COMP(j).Fs4tile = nansum(tmpE,3);
        tmpF = cat(3,COMP(j).Fs5tile,Fs5_gen); COMP(j).Fs5tile = nansum(tmpF,3);        
        tmpH = cat(3,COMP(j).Fs6tile,Fs6_gen); COMP(j).Fs6tile = nansum(tmpH,3);
        tmpI = cat(3,COMP(j).Fs7tile,Fs7_gen); COMP(j).Fs7tile = nansum(tmpI,3);
        tmpJ = cat(3,COMP(j).Fs8tile,Fs8_gen); COMP(j).Fs8tile = nansum(tmpJ,3);
        tmpG = cat(3,COMP(j).Zvegtile,Zveg_gen); COMP(j).Zvegtile = nansum(tmpG,3);
        
        % aggregate the lumped flux
        tmpA = cat(3,COMP(j).Fs0tile_abs,CH4_genL); COMP(j).Fs0tile_abs = nansum(tmpA,3);
        tmpB = cat(3,COMP(j).Fs1tile_abs,s1_genL); COMP(j).Fs1tile_abs = nansum(tmpB,3);
        tmpC = cat(3,COMP(j).Fs2tile_abs,s2_genL); COMP(j).Fs2tile_abs = nansum(tmpC,3);
        tmpD = cat(3,COMP(j).Fs3tile_abs,s3_genL); COMP(j).Fs3tile_abs = nansum(tmpD,3);

        %% Calculate weights
        
        % Calculate the weight of each cell
        tmpW = cat(3,COMP(j).footWeight,(FPall(j).fCrop)); COMP(j).footWeight = nansum(tmpW,3);
        % Calculate the number of times each tile is counted
        tmpS = cat(3,COMP(j).footCount,~isnan(FPall(j).fCrop)); COMP(j).footCount = nansum(tmpS,3);
        
        % Do the same but for those time where FCH4 was available
        if ~isnan(FCH4(i))
        tmpW = cat(3,COMP(j).footWeightS0,(FPall(j).fCrop)); COMP(j).footWeightS0 = nansum(tmpW,3);    
        tmpS = cat(3,COMP(j).footCountS0,~isnan(FPall(j).fCrop)); COMP(j).footCountS0 = nansum(tmpS,3);
        end
        
        end
    end  
        

else
     disp (['Did not run footprint for day ' num2str(DOY(i)) ', time: ' num2str(endH(i))])
end

end


end