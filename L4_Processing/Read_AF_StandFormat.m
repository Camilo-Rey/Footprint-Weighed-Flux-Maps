
% Read AMF data
% Function to transform csv files from Ameriflux in Standard Format to
% Matlab files
% Created by Camilo Rey-Sanchez.

function [AF]=Read_AF_StandFormat(filename,filepath)


opts = detectImportOptions(strcat(filepath,filename));
Data = readtable(strcat(filepath,filename),opts);

% Mdate
datestring = num2str(Data.TIMESTAMP_START);
year = str2num(datestring(:,1:4));month = str2num(datestring(:,5:6));
day = str2num(datestring(:,7:8));hr = str2num(datestring(:,9:10));min = str2num(datestring(:,11:12));
Mdate = datenum(year, month, day, hr, min, 0);
dateplot=datetime(datevec(Mdate));

% Index by year
ind=find(~isnan(year));
AF.Mdate=Mdate(ind);
AF.dateplot=dateplot(ind);
dimNames = Data.Properties.VariableNames;

for j=3:length(dimNames)

name=dimNames(j);
AF.(name{1}) = table2array(Data(ind,j));
AF.(name{1})(AF.(name{1})==-9999)=nan; % Remove -9999

end

%figure;plot(AF.dateplot,AF.P);title('Prec');
%save(['AF_' site],'AF')