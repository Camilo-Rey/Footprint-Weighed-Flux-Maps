%% Estimate h based on Pennypacker and Baldocchi (2015) paper using values obtained above
function [h_dsk, zh_final]=CalculateAerodynamicCanopyHeight(Mdate,ustar,z,L,u,daysAVG,plotYN)
k = 0.4; 
d_h = 0.6; %Default from paper
z0_h = 0.1; %Default from paper
Dplot=datetime(datevec(Mdate));
% Define critaria for calculating h

% Most common filter
ind = find(ustar > 0.2 & ustar < 0.3 & abs(z./L) < 0.01 & u>1.5);

% Short section
% ind = find(ustar > 0.15 & ustar < 0.4 & abs(z./L) < 0.05 & u>1.5);
% disp('WARNING, using short data set filters for aerodynamic canopy height')

%ind = find(ustar > 0.25 & ustar < 0.35 & abs(z./L) < 0.01 & (WD > 270 & WD < 315));

h = NaN(size(ustar));
h(ind) = z(ind)./(d_h+z0_h.*exp(k.*u(ind)./ustar(ind)));


% Calculate daily mean height using smoothing and despike one last time
Days=day(Dplot,'dayofyear');
[h_despike,~] = despike(h(ind),7,0.5);
if length(unique(Days))>30 % higher than 30 days of data
%plot([h(ind),h_despike],'.')
navg=daysAVG;% number of windows to use for average
hsmooth_dsk=smooth(h_despike,navg);
h_dsk = NaN(size(ustar));
h_dsk(ind) = hsmooth_dsk;
% Interpolate daily measurements back to 30 min values
zh_final = interp1(Mdate(~isnan(h_dsk)),h_dsk(~isnan(h_dsk)),Mdate,'linear','extrap');
% Fix the interpolation at the end
fff=find(~isnan(h_dsk));zh_final(1:fff(1))=zh_final(fff(1));
% Fix the interpolation at the end
zh_final(fff(end):end)=zh_final(fff(end));
else
disp('Too little data, z_veg will be equal to the average')
zh_final=ones(length(h),1)*nanmean(h_despike);
end

if plotYN==1 && length(unique(Days))>4 % higher than 4 days of data
figure(100);
plot(Dplot,h,'.k','MarkerSize',14)
hold on
plot(Dplot, h_dsk,'.b','MarkerSize',20);
plot(Dplot, zh_final,'-b');
legend('initial','Despiked','smoothed')
ylabel('Aerodynamic canopy height (m)')
title(['Smoothing based on ' num2str(daysAVG) '-day averages of canopy height'])
end
end


