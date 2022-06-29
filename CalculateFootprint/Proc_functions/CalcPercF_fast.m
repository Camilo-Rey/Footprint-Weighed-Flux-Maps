function[tempMap2]=CalcPercF_fast(footD)
    percOpts  = sort(unique(footD),'descend');
    tempMap2=nan(size(footD));
        for i=1:length(percOpts)
        curPerc=percOpts(i);
        good=footD>=curPerc;
        score=sum(sum(footD(good)));
        tempMap2(good & isnan(tempMap2))=score;
        if score>0.99*nansum(nansum(footD)) % If the score is higher than 99% of the total, break the process to save resources
        break
        end    
        end
end