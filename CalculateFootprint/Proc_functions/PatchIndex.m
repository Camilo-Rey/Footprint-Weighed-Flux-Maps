function [PI] = PatchIndex(patch)

PI=unique(patch);

if length(PI)==2 && min(PI)==0 && max(PI)==1
    PI=1;
end

end
    