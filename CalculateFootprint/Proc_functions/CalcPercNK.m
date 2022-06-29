function [FFP]= CalcPercNK(footC,FX2,FY2,Rl,contV);
 f_array = reshape(footC,1,size(footC,1)*size(footC,2));
       f_sort  = sort(f_array,'descend');
       f_sort  = f_sort(~isnan(f_sort));
       f_cum   = cumsum(f_sort);
       
       GG=nansum(f_sort);
       if GG<0.9;contV(2)=90*GG;end
       for i = 1:length(contV)
           f_diff     = abs(f_cum-contV(i)/100);
           [~, ind_r] = min(f_diff);
           fr         = f_sort(ind_r);
           contour_r  = contourc(FX2(1,:),FY2(Rl:-1:1,1),footC,[fr fr]);    
           % Contourc adds info on level and number of vertices - replace with NaN
           ind_nan = contour_r(1,:)==fr;
           contour_r(:,ind_nan) = NaN;
           % Decrease number of digits and sort/unique
               contour_r = round(contour_r,1);
                           
               % Fill output structure               
               FFP(i).fr  = fr;
               FFP(i).xr  = contour_r(1,:);
               FFP(i).yr  = contour_r(2,:);
       end
end
