%% Code to create Footprint-weighed Flux maps
%  This code takes the input from CalculateFootprint and plots a flux map
%  of the desired scalar.
% Created by Camilo Rey-Sanchez
% Last modified: Jun 30 2022

site='NC4';

DataDir=(['Q:\My Drive\NC-State\Ameriflux_data\FW_maps\' site '\Footprint_Output\continuous\']);% Modify the data directory where the footprint output was saved

FileToLoad= 'NC4-KM-2017274-0000-to-2017288-2400-day-50-80test';% Modify the name of the file accordingly    

load([DataDir FileToLoad])    

FPall(1).name='Hsieh';
FPall(2).name='Kljun';
FPall(3).name='K&M';


%% Surface Flux Map

cutFP=80;

figure('position',[50,100,1100,370])
 for j=3 % Number 3 is K&M model     

    % CH4 flux 
    subplot(1,3,1)
    Lflux=COMP(j).Fs0tile./COMP(j).footWeightS0;
    Lflux(FG(j).percF>cutFP)=nan; Lflux(isnan(FG(j).percF))=nan;   
    pcolor(FX2,FY2(RT:-1:1,:),Lflux);
    shading('flat'); colorbar('Location','south');  
    colormap(jet)     
    xlabel('x [m]');  ylabel('y [m]')
    title('CH_4 flux (nmol m^{-2} s^{-1})')    
    set(gca,'fontsize',12)
 
    
    % CO2 flux 
    subplot(1,3,2)
    Lflux=COMP(j).Fs1tile./COMP(j).footWeight;
    Lflux(FG(j).percF>cutFP)=nan; Lflux(isnan(FG(j).percF))=nan; 
    pcolor(FX2,FY2(RT:-1:1,:),Lflux);
    shading('flat'); colorbar('Location','south');     
    xlabel('x [m]'); 
    title(opt.FS1) 
    set(gca,'fontsize',12)
    
     % H2O flux 
    subplot(1,3,3)
    Lflux=COMP(j).Fs2tile./COMP(j).footWeight;
    Lflux(FG(j).percF>cutFP)=nan; Lflux(isnan(FG(j).percF))=nan;   
    pcolor(FX2,FY2(RT:-1:1,:),Lflux);
    shading('flat'); colorbar('Location','south');     
    xlabel('x [m]');   
    title(opt.FS2)   
    set(gca,'fontsize',12)     

    
   end
 
 
