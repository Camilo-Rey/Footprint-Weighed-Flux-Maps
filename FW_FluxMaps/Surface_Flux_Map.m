%% Code to create Footprint-weighed Flux maps
%  This code takes the input from CalculateFootprint and plots a flux map
%  of the desired scalar.
% Created by Camilo Rey-Sanchez
% Last modified: Jun 30 2022

site='NC4';

PLT=6;% Odd are short frames, even are large frames

DataDir=(['Q:\My Drive\UC-Berkeley\Footprint_share\Footprint_Output\' site '\continuous\']);% Modify the data directory where the footprint output was saved

FileToLoad= 'NC4-KM-2019168-0000-to-2019192-2400-day-50-80';% Modify the name of the file accordingly    

load([DataDir FileToLoad])    

FPall(1).name='Hsieh';
FPall(2).name='Kljun';
FPall(3).name='K&M';


%% Surface Flux Map

cutFP=80;

figure('position',[50,100,1100,370])
 for j=3
     

    % CH4 flux 
    subplot(1,3,1)
    Lflux=COMP(j).CH4tile./COMP(j).footWeightCH4;
    Lflux(FG(j).percF>cutFP)=nan; Lflux(isnan(FG(j).percF))=nan;   
    pcolor(FX2,FY2(RT:-1:1,:),Lflux);
    shading('flat'); colorbar('Location','south');  
    colormap(jet)     
    xlabel('x [m]');  
    title('CH_4 flux (nmol m^{-2} s^{-1})')    
    set(gca,'fontsize',12)
 
    
    % CO2 flux 
    subplot(1,3,2)
    Lflux=COMP(j).CO2tile./COMP(j).footWeight;
    Lflux(FG(j).percF>cutFP)=nan; Lflux(isnan(FG(j).percF))=nan; 
    pcolor(FX2,FY2(RT:-1:1,:),Lflux);
    shading('flat'); colorbar('Location','south');     
    xlabel('x [m]'); 
    title('CO_2 flux (\mumol m^{-2} s^{-1})') 
    set(gca,'fontsize',12)
    
     % H2O flux 
    subplot(1,3,3)
    Lflux=COMP(j).H2Otile./COMP(j).footWeight;
    Lflux(FG(j).percF>cutFP)=nan; Lflux(isnan(FG(j).percF))=nan;   
    pcolor(FX2,FY2(RT:-1:1,:),Lflux);
    shading('flat'); colorbar('Location','south');     
    xlabel('x [m]');   
    title('H_2O flux (mmol m^{-2} s^{-1})')   
    set(gca,'fontsize',12)     

    
   % sgtitle([FPall(j).name ' DOY ' opt.start(9:11) 'to' opt.end2(9:11) ' Year ' opt.start(1:4)])  
 end
 
 
