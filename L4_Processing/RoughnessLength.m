function [d,z0,n,stat]=RoughnessLength(TowerHeight,CanopyHeight,OL,ustar,wind_speed,daynight)
%RoughnessLength: Calculates displacement height and roughness length from
%the parameters provided
% Based on Mauer et al (2013) from Gil Bohrer's lab.

lim=-1;

z=TowerHeight;
hei=CanopyHeight;

good = ((z./OL>lim & z./OL<0) & ~isnan(wind_speed) & ~isnan(ustar) );

z_star=hei*2;

USTAR=ustar(good);
L=OL(good);
ubar=wind_speed(good);

if sum(good)>20% previosuly 200
    
    
    s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR','MaxFunEvals',1000,'MaxIter',1000,'Display','Off','Upper',[hei hei],'Lower',[0 0],'StartPoint',[5 1.5]);

    f1 = fittype('USTAR/0.4*(log((z-a)/b)-(2*log((1+(1-16*(z-a)/x)^0.25)/2)+log((1+(1-16*(z-a)/x)^0.5)/2)-2*atan((1-16*(z-a)/x)^0.25)+pi/2)+(2*log((1+(1-16*(b)/x)^0.25)/2)+log((1+(1-16*(b)/x)^0.5)/2)-2*atan((1-16*(b)/x)^0.25)+pi/2)+((1-16*(z-a)/x)^(-1/4))*((1+0.5/2.59/((z-a)/(z_star-a)))*(z-a)/x)/1.5*log(1+1.5/(2.59*((z-a)/(z_star-a))))*exp(-2.59*((z-a)/(z_star-a))))',...           
        'problem',{'z','USTAR','z_star'},'options',s);

    [P,gof] = fit(L,ubar,f1,'problem',{z,USTAR,z_star})

    d=P.a; %Displacement height
    z0=P.b; %Roughness length
    n=sum(good);
    stat=gof;
else
    d=nan;
    z0=nan;
    n=nan;
    stat=nan;
end
