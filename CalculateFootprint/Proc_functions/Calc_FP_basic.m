% Footprint calculation using multiple models
% Created by Camilo Rey-Sanchez. (A big chunk of this code is legacy from
% my work in the lab of Gil Bohrer, so acknowledgement to him and lab
% members, particularly Tim Morin are due here.
%
% Inputs:
%        ustar   = friction velocity [m/s]
%
%        Lo      = Obukhov stability length [m]
%            
%        sv      = fluctuation of lateral wind [m/s] (approx. 2 times ustar)
%                       
%        zo      = momentum roughness length [m]
%
%        zH      = height of the measurement [m]
%          
%        d - displacement height
%
%       windir - wind direction, from north, degrees. Where wind is coming
%       FROM, not going to
%       
%        FX - meshgrid, x coordinates, of mask map
%
%        FY - meshgrid, y coordinates, of mask map
%
%        F - mask map on FX,FY coordinates

function [FPall,pxR]=...
    Calc_FP_basic(ustar,Lo,sv,zo,zH, d, windir, plotYN, FX, FY, PatchMap,Cmap,...
    h,contourMarks,Tdist,Twindir,ubar,models,mscale)

%% Overall Parameters

zm=zH-d; % Effective heigth, tower heigth minus displacement
Lx=100*zm;% The initial length of the along-wind footprint is equal to 150the heigth of the measurement
k=0.4; % von Karman constant
zL=zm/Lo; % Stability coefficieny
[Txmax,Tymax]=rotateToWind(Tdist,0,Twindir);
[FX2, FY2] = meshgrid(FX,FY);% The coordinates of our map
[Fxprime,Fyprime] = rotateToWind(FX2,FY2,windir);%     
Rl= length(FY);

% Plotting footprint parameters
    bin=0.5;
    ywidth0=floor(zo*(0.3*(sv/ustar).*(Lx./zo).^0.86)/1.5);% Defining the width of the lateral dispersion
    ywidth=ywidth0*3;% 3 times the width 

% graph  
if plotYN==1
    figure;
    if ~isnan(PatchMap)
    pcolor(FX2,FY2(Rl:-1:1,:),PatchMap*100);shading('flat');% This is the same as using flidud, which flips the rows so that the first one is the last one
    colormap(flipud(parula))
    colorbar
    else
    pcolor(FX2,FY2(Rl:-1:1,:),nan(length(FX),length(FY)));shading('flat');% This is the same as using flidud, which flips the rows so that the first one is the last one
    end    
    hold on; xlabel('x [m]');  ylabel('y [m]')
end

FPall(1).name='Hsieh';
FPall(2).name='Kljun';
FPall(3).name='K&M';

zu=zm*(log(zm/zo)-1+zo/zm);% New heigth scale. Eq. 13.5 Hsieh (2000)   

% Coefficients D and P according to stability
    P=[0.59 1 1.33];                  
    D=[0.28 0.97 2.44];

    stab=zu/Lo;
    thresh=0.04;

    if stab<-thresh
        ii=1;
    elseif abs(stab)<thresh
        ii=2;
    elseif stab>thresh
        ii=3;
    end
    D1=D(ii);
    P1=P(ii);
   

    F2H=(D1/0.105/k^2)*(zm^(-1)*abs(Lo)^(1-P1)*zu^(P1)); % Fetch to height ratio. Eq 20, Hsieh (2000)
    Xm=ceil(F2H*zm);
%% Hsieh model

if models(1)==1
tic    
    if mscale==0 % linear
        bin=max(1,round(F2H*zm/500,0)); % length of the bin in meters. 1 m is the minimum
        x=eps:bin:Xm; x(x==0)=0.001;  % Sequence of point along the main axis

    elseif mscale==1 % exponential
        bin=0.1;
        x=exp(0:bin:log(Xm*1.1));% 10 % higher than F2H to encompass a bit more area
    end

    T2=(-1/k/k)*(D1*zu^P1*abs(Lo)^(1-P1))./(x);  % Eq 17. Hsieh (2000)   
    Fp=-(T2./x).*exp(T2);% Cross-wind integrated (along wind) footprint function

% Detto (2006) 2D expansion

    nn=length(x);% The number of bin increments in the distance Xm
    sy=zo*0.3*(sv/ustar).*(x./zo).^0.86;% Eq B4 Detto. Standard deviation of the cross-wind (besed on the standard deviation of lateral fluctuation sv)

    
    if mscale==0 % linear
        y=(-ywidth:bin:ywidth);   
    elseif mscale==1 % exponential
    	y=[flip(-exp(0:bin:log(ywidth))),0, exp(0:bin:log(ywidth))];
    end
    
    footo=nan(nn,length(y)); % Final source function matrix with nn rows (along-wind) and y columns (cross-wind)
    for i=1:nn
        footo(i,:)=Fp(i)*normpdf(y,0,sy(i));% This normpdf function is equivalent to the equation B3 in Detto (2006)
    end

% Unrotated footprint
    foot=footo';
    % The sum of foot must be close to unity after integrating using the correct pixel
    % size (usually 0.5 * 0.5 m = 0.25 m2)

%% Arranging and rotating matrices

if ~all(y==0)
    [xx,yy] = meshgrid(x,y);% The coordinates of foot after transposing
          
    footRotate = griddata(xx,yy,foot,Fxprime,Fyprime); % This puts f(x,y) in the right direction and inside the map grid

    % binnorm first creates a factor to adjust for pixel size (which is usally
    % not equal to one). Secondly it creates a factor to normalized the rotated
    % footprint so that it is equal to the sum of the original unrotated
    % footprint (which must be close to unity)
    %Still not sure why the rotation causes such differences.
    
    if mscale==0
    dx=(x(2)-x(1));
    dy=(y(2)-y(1));
    binnorm = (nansum(nansum(foot))*dx*dy)/nansum(nansum(footRotate));% The size of the bin needs to be normalized 
    footHsieh=footRotate.*binnorm;% This is the corrected footprint that takes into account the change in scale

    elseif mscale==1
    dx=[x(1),(x(2:end)-x(1:end-1))];
    dy=[abs(y(2)-y(1)),abs(y(2:end)-y(1:end-1))];
    [dxx,dyy] = meshgrid(dx,dy);
    DD=foot.*dxx.*dyy;
    binnorm=nansum(nansum(DD))/nansum(nansum(footRotate));
    footHsieh=footRotate.*binnorm;% This is the corrected footprint that takes into account the change in scale
    nansum(nansum(footHsieh)); % This should be close to 0.8    
    end
    
    
%% Finding other points of interest

% Looking at the distance with the maximum contribution:
    HXmax=D1*zu^P1*abs(Lo)^(1-P1)/(2*0.4^2); % Eq 19 (Hsieh, 2000)
    [~,hh]=max(Fp);
    HXmax2=x(hh);% another alternative to find the maximum
    [Hxmax,Hymax]=rotateToWind(HXmax,0,windir);

%% Calculate Percent Footprint and plot  
 
    [FFPo]= CalcPercNK(footHsieh,FX2,FY2,Rl,contourMarks);

% PLOT

if plotYN==1
   figure(gcf);%   
     for i=1:length(contourMarks)          
           LG(1).p= plot(FFPo(i).xr,FFPo(i).yr,'k');hold on;
     end  
     LG(1).o=plot(Hxmax,Hymax,'ok');
end

fCropH=footHsieh;
fCropH(fCropH<FFPo(end).fr)=nan;

FPall(1).foot=footHsieh;
FPall(1).fCrop=fCropH;
FPall(1).Xmax=HXmax;

end
toc
end

%--------------------------------------------------------------------------
%% Kljun model
%--------------------------------------------------------------------------

if models(2)==1
tic    
if zL>=-15.5

[FFP,~]=calc_footprint_FFP_2(zm,zo,ubar,h,Lo,sv,ustar,...
    'wind_dir',windir,'fig',0,'r',contourMarks,'nx',600,...
    'rslayer',1);

    
%% Arranging matrices

    footRotateK = griddata(FFP(1).x_2d,FFP(1).y_2d,FFP(1).f_2d,FX2,FY2(Rl:-1:1,:)); % This puts f(x,y) in the right direction and inside the map grid
    
    xpp=FFP(1).x_ci;ypp=FFP(1).y;
    binnorm2 = (nansum(nansum(FFP(1).f_2d))*(xpp(2)-xpp(1))*(ypp(2)-ypp(1)))...
        /nansum(nansum(footRotateK));% The size of the bin needs to be normalized 
    
    footKljun=footRotateK.*binnorm2;%
    Ksum=nansum(nansum(footKljun));
    [Kxmax,Kymax]=rotateToWind(FFP(1).x_ci_max,0,windir);
    
    if plotYN==1
    figure(gcf);%
     for i=1:length(contourMarks)          
        LG(2).p= plot(FFP(i).xr,FFP(i).yr,'r');hold on;
     end  
     hold on   
     LG(2).o=plot(Kxmax,Kymax,'or'); 
    end

fCropK=footKljun;
fCropK(fCropK<FFP(end).fr)=nan;

FPall(2).foot=footKljun;
FPall(2).fCrop=fCropK;
FPall(2).Xmax=FFP(1).x_ci_max;
FPall(2).f80=FFP.x_80;


else    
    FPall(2).Xmax=nan;
    FPall(2).foot=nan(length(FX),length(FY));
    FPall(2).fCrop=nan(length(FX),length(FY));
    footKljun=nan(length(FX),length(FY));
end
toc
end
%--------------------------------------------------------------------------
%% Korman and Meixner
%--------------------------------------------------------------------------

if models(3)==1
tic    
% if mscale==0% linear    
%     x=eps:bin:Xm; x(x==0)=0.001; xKM=x; % re-initialize the length
% elseif mscale==1;    
%     x=exp(0:bin:log(Xm*1.1));% 10 % higher than F2H to encompass a bit more area
% end

    if mscale==0 % linear
        bin=max(1,round(F2H*zm/500,0)); % length of the bin in meters. 1 m is the minimum
        x=eps:bin:Xm; x(x==0)=0.001;  % Sequence of point along the main axis

    elseif mscale==1 % exponential
        bin=0.1;
        x=exp(0:bin:log(Xm*1.1));% 10 % higher than F2H to encompass a bit more area
    end

 nn=length(x);% The number of bin increments in the distance Xm

% Parameters
    if zL>0;
        zt=0;
        phim=1+5*zL;% Stability coefficient for wind shear
        phic=1+5*zL;% Stability coefficient for eddy diffusivity of any scalar(k). Assumed to equal to diffusivity for heat (van Ulden,1978)
        psim=5 * zL;
        n=1/phic;
    else
        zt=(1 - 16 * zL)^0.25;
        phim=(1-16*zL)^-0.25;% Stability coefficient for wind shear
        phic=(1-16*zL)^-0.5;% Stability coefficient for eddy diffusivity of any scalar(k). Assumed to equal to diffusivity for heat (van Ulden,1978)
        psim=-2 * log((1 + zt)/2) - log((1 + zt^2)/2) + 2 * atan(zt) - pi/2;
        n=(1 - 24 * zL)/(1 - 16 * zL);
    end
    uz=ustar/k*(log(zm/zo)+psim);
    eddydif = k * ustar * zm/phic;% use the stability coefficient for eddy diffusivity of any scalar
    
    m = ustar * phim/(0.41 * uz);
    r=2+m-n;% Shape factor (Van Ulden, 1978)
    mu=(1+m)/r;
    fgamma=gamma(mu);
    alpu = uz/(zm^m);% Eq 11 K&M 
    alpk = eddydif/(zm^n);% Eq 11 K&M 
    xi= alpu*zm^r/(r^2*alpk);% Eq 19 K&M 
    FpKM=(1/fgamma) .* (xi^mu)./(x.^(1 + mu)) .* exp(-xi./x);% Eq 21 K&M 

% lateral dispersion

    if mscale==0 % linear
        y=(-ywidth:bin:ywidth);% two times the width    
    elseif mscale==1 % exponential
    	y=[flip(-exp(0:bin:log(ywidth))),0, exp(0:bin:log(ywidth))];
    end
   
     [xx,yy] = meshgrid(x,y);% The coordinates of foot after transposing
     
    footKMo=nan(nn,length(y)); % Final source function matrix with nn rows (along-wind) and y columns (cross-wind)
    
    %footKMold=nan(nn,length(y)); % Final source function matrix with nn rows (along-wind) and y columns (cross-wind)
%     for i=1:nn
%         footKMold(i,:)=FpKM(i)*normpdf(y,0,sy(i));% This normpdf function is equivalent to the equation B3 in Detto (2006)
%     end

% lateral dispersion new

    uPlume = (gamma(mu)/gamma(1/r));
    uPlume = uPlume * ((r^2 * alpk/alpu)^(m/r)); % Eq 18. K&M (2000)
    auPlume = uPlume * (alpu * (x.^(m/r)));
    sigmaY = sv * x./auPlume;
    
    for i=1:nn
    footKMo(i,:)=FpKM(i)*normpdf(y,0,sigmaY(i));% This normpdf function is equivalent to the equation B3 in Detto (2006)
    end

%     footprint = Pfp/sum(sum(Pfp))
%     returnList = list(footprint = footprint, FPe = FPe, FPn = FPn)
 
%% Arranging and rotating matrices

% Unrotated footprint
    footF=footKMo';

if ~all(y==0)
    footRotateKM = griddata(xx,yy,footF,Fxprime,Fyprime); % This puts f(x,y) in the right direction and inside the map grid

    
    if mscale==0
    dx=(x(2)-x(1));
    dy=(y(2)-y(1));
    binnorm = (nansum(nansum(footF))*dx*dy)/nansum(nansum(footRotateKM));% The size of the bin needs to be normalized 
    footKM=footRotateKM.*binnorm;% This is the corrected footprint that takes into account the change in scale
    elseif mscale==1
    dx=[x(1),(x(2:end)-x(1:end-1))];
    dy=[abs(y(2)-y(1)),abs(y(2:end)-y(1:end-1))];
    [dxx,dyy] = meshgrid(dx,dy);
    DD=footF.*dxx.*dyy;
    binnorm=nansum(nansum(DD))/nansum(nansum(footRotateKM));
    footKM=footRotateKM.*binnorm;% This is the corrected footprint that takes into account the change in scale
    nansum(nansum(footKM)); % This should be close to 0.8    
    end

% Distance of the maximum footprint contribution Xmax
    KMXmax=xi/(1+mu);
    [KMxmax,KMymax]=rotateToWind(KMXmax,0,windir);

% calculate coordinates for countours   
    [FFP3]= CalcPercNK(footKM,FX2,FY2,Rl,contourMarks);

% PLOT

if plotYN==1
        for i=1:length(contourMarks)          
        LG(3).p= plot(FFP3(i).xr,FFP3(i).yr,'b');hold on;        
        end
   hold on   
   LG(3).p=plot(KMxmax,KMymax,'ob');  
   
end

fCropKM=footKM;
fCropKM(fCropKM<FFP3(end).fr)=nan;

FPall(3).foot=footKM;
FPall(3).fCrop=fCropKM;
FPall(3).Xmax=KMXmax;

end
toc
end

if plotYN==1
   to=plot(Txmax,Tymax,'ob');   
   LH=[];
   LLAB={'Legend'};
   for i=1:3
       if models(i)==1
       LH=[LH,LG(i).p];
       LLAB=[LLAB,[FPall(i).name ' ' num2str(contourMarks(end)) '%']];
       LLAB=[LLAB,[FPall(i).name ' Max Source']];
       end
   end
   LLAB=[LLAB,'Point of interest'];
   legend(LLAB,'Location','Best')
end

%% Calculate the cell in the footprint from a given distance (The release in this case)
           
    [~,closestIndexX] = min(abs(FX2(1,:)-Txmax));
    [~,closestIndexY] = min(abs(FY2(Rl:-1:1,1)-Tymax));              

        if models(1)==1;pxR.H=footHsieh(closestIndexY,closestIndexX);end
        if models(2)==1;pxR.K=footKljun(closestIndexY,closestIndexX);end
        if models(3)==1;pxR.KM=footKM(closestIndexY,closestIndexX);end

%          center=[find(FX2(1,:)==0),find(FY2(:,1)==0)];
%             nCol=round((Txmax/pxSize));
%     nRow=round(Tymax/pxSize);
%     cM2=center+[-nRow,nCol];
%     if cM2(1)>0 && cM2(2)>0
%         
%         if models(1)==1;pxR.H=footHsieh(cM2(1),cM2(2));end
%         if models(2)==1;pxR.K=footKljun(cM2(1),cM2(2));end
%         if models(3)==1;pxR.KM=footKM(cM2(1),cM2(2));end
%     else
%     pxR.H=nan;   pxR.K=nan;   pxR.KM=nan;
%     end
    
    
    
% Plot one dimensional
if plotYN==1
%     figure(2); 
%     if models(1)==1;plot(x,Fp);hold on;end
%     if models(2)==1;plot(xKljun,kfy);hold on;end
%     if models(3)==1;plot(xKM,FpKM);hold on;end
%     xlabel('Distance from tower, x (m)')
%     ylabel('Cross-wind integrated footprint')
%     legend('Hsieh 2000','Kljun 2015','Korman & Meixner 2001')
%     hold off;

end

%% Filling with nans
for yt=1:3
if models(yt)~=1
    FPall(yt).Xmax=nan;
    FPall(yt).foot=nan(length(FX),length(FY));
    FPall(yt).fCrop=nan(length(FX),length(FY));
end
end
%% This part adds up all the values of the source function that come from
% different i land cover types

    
PI = PatchIndex(PatchMap); % If there is a land cover map, this will be used.
if length(PI)<20 % More than 20 means it is probably not a land cover map
for i=1:length(PI)
    if models(1)==1;FPall(1).PatchSum(i) = nansum(nansum(footHsieh.*(PatchMap==PI(i))));end% By multiplying by the size of the bin in m2, we are normalizing by area, but this is not necessary
    if models(2)==1;FPall(2).PatchSum(i) = nansum(nansum(footKljun.*(PatchMap==PI(i))));end% By multiplying by the size of the bin in m2, we are normalizing by area, but this is not necessary
    if models(3)==1;FPall(3).PatchSum(i) = nansum(nansum(footKM.*(PatchMap==PI(i))));end% By multiplying by the size of the bin in m2, we are normalizing by area, 
end
end


end




