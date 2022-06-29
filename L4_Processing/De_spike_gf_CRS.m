
function [wc2,wc3] = De_spike_gf_CRS(wc,m,tr,npts_day,min,max)
%% Function to 1) Despike a dataset based on a triming value (tr) and
% 2) Gapfill based on a linear interpolation from a 3D grid data.
% The first dimension is the diurnal cycle, the second dimension is a vector of daily trends.
% The third dimension (Z) is the dataset
% Created by Gil Bohrer lab. Commented and Improved by Camilo Rey on Jan 2, 2018.

%% Inputs:
% wc: dataset
% m: Number of days of the time window to evaluate the spike
% tr: trimming value. How many std from the mean of ma days is considered a spike
% npts_day: Number of points in a day
%% Outputs:
% wc2: Despiked vector
% wc3: Gap-filled vector

if nargin<5
   min=min(wc);
   max=max(wc);
end

out=wc>max | wc<min;
wc(out)=nan;

NN=length(wc);                         % Length of vector to despike
nw=npts_day*m ;                      % Number of points in the time window
R=floor(NN/nw) ;                 % Number of realizations for the despike

WC=wc; %Original vector
wc2=wc; % Despiked vector

%% DESPIKING

 for y=1:R
        Touse=[(y-1)*nw+1:y*nw]';
        % The last realization may not be of the same lenght of nw so we
        % grab wahtever number of points until the end of the dataset:
        if Touse(end)>NN-nw
           Touse=[(y-1)*nw+1:NN]';
        end
 
 hh=find(~isnan(wc(Touse)));% Find the non-nan values

 if length(hh)>0.1*nw
 a=polyfit(Touse(hh),wc(Touse(hh)),1);%find a linear trend in the length of the dataset (e.g. seasonal trend)
 wcs=wc(Touse)-polyval(a,Touse);% remove the seasonal trend from dataset
% figure;plot(wc(Touse));hold on;plot(polyval(a,Touse));
% title(['linear trend in window ' num2str(y) ' out of ' num2str(R) ' realizations'])
 
 % Standard is a nested function at the end of this script. It calculates
 % the median for each half hour and calculates the difference from each
 % half hour to that median
 %[wcsB,~]=standard(wcs,npts_day);% stnd is a vector with the median. wcs is the difference 
 wcsC=red(wcs); %red is a function that normalizes all the points to mean o and std 1
  %figure;plot(wc);hold on;plot(wcs);plot(wcsB);plot(wcsC);
  %legend('raw data','de-trended','standardized to median','standardized to mean 0')
 
 %These are the differences between adjacent points
 dwcs1=[diff(wcsC); 0];
 dwcs2=[0; diff(wcsC)];
 %There are 3 conditionals to find spikes: higher than tr standard
 %deviations from point above. 2) higher than tr std from point below and
 %3) higher than tr stds from the median
 bad=find(abs(red(dwcs1))>tr | abs(red(dwcs2))>tr |abs(wcsC)>tr | wc(Touse)>max | wc(Touse)<min);
     nbad=length(bad);
     disp([num2str(nbad) 'spikes removed out of ' num2str(NN)])
     wc2(Touse(bad))=NaN;
 %figure;plot(wc2(Touse));
 else
 disp('Too many nans to perform a despike or gapfill procedure')
nbad=48; 
wc2(Touse)=WC(Touse);
wc3(Touse)=WC(Touse);
 end
 end


%% GAPFILLING
% Gap-filling with 2-D (griddata) look-up table (along time and day)

% SOmetimes the first value of a day is nan:
if length(hh)>0.1*nw

if nbad==0
wc3=wc2;
else
yy=reshape(wc2,npts_day,NN/npts_day);
for j=1:size(yy,2)
if isnan(yy(1,j))==1 & size(yy,1)>1
    yy(1,j)=yy(2,j);
end
end
yyy=reshape(yy,NN,1);
run=find(~isnan(yyy));
gap1=find(isnan(yyy));

nday=NN/npts_day;
XI(1,:)=1:nday;
YI(:,1)=1:npts_day;
X=repmat(XI,npts_day,1);

%disp(['X:' num2str(size(X,1)) ' Y:' num2str(size(X,2))
X=reshape(X,1,nday*npts_day);
Y=repmat(YI,nday,1);
Z=yyy;

method='linear';
ZI = griddata(X(run),Y(run),Z(run),XI,YI,method);
wc3=reshape(ZI,NN,1);
wc3(run)=wc2(run);

wc3=wc2;
%extrapolate using nearest
run=find(~isnan(wc3));
gap=find(isnan(wc3));
if ~isempty(gap)
Z=wc3;
method='nearest';
ZI = griddata(X(run),Y(run),Z(run),XI,YI,method);
Z=reshape(ZI,NN,1);
wc3(gap)=Z(gap);

end

figure;
plot(wc3,'k');hold on
plot(wc2)
plot(gap,WC(gap),'k+')
legend('gap filled','despiked','gapfilled 2')
title([num2str(length(gap1)) 'gaps filled. ' num2str(sum(isnan(wc3))) ' nans remaining in vector'])
disp([num2str(length(gap1)) 'gaps filled. ' num2str(sum(isnan(wc3))) ' nans remaining in vector'])
end
else
wc2=WC;
wc3=WC;

end
%% Nested functions

%standarizzation of time series
% function [ys,stnd] = standard(y,period)
% 
% n=length(y);
% y1(:,1)=y;
% 
% Yv=reshape(y1,period,n/period);
% mY=nanmedian(Yv,2);% The median for that half hour of the day.
% 
% stnd=repmat(mY,n/period,1); 
% ys=(y-stnd);% Difference from the median
% 
% end

% use=find(isnan(ys)==1);
% ys(use)=0;
%end

% We dont use this function, but I'll be it here in case we need it.
%find gaps<=M and fill with spline iterpolant

function y_gf=gap_filler(x,y,M)

if isnan(y(1))==1;
    y(1)=y(2);
end

y_gf=y;
use=find(isnan(y)==0);
gap_M=[];

for i=1:length(y)-M-1
    if isnan(y(i))==1 & sum(isnan(y(i:i+M+1)))<M & isnan(y(i-1))==0
       b=find(isnan(y(i:end))==0);
       gap_M=[gap_M i:i+b(1)-2];
       i=i+b(1)-1;
    end
end
  

y_gf(gap_M)=interp1(x(use),y(use),x(gap_M),'spline');

end

%normilize the variable x with zero mean and unit vatiance

function xn=red(x)


     if nanstd(x)==0
       xn=x*0;
     else
       xn=(x-nanmean(x))/nanstd(x);
     end
end

end
