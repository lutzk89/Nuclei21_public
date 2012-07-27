function posP = invPosPeriodicLE( pos, rim1, dims )

if(nargin <= 2)
    dims = 1:size(pos,2);
end
rim2 = 1-rim1;

posP = pos;
for nu=dims %cycle through all boundaries
    posAdd1 = posP(posP(:,nu)<rim1(nu),:); 
    posAdd2 = posP(posP(:,nu)>rim2(nu),:); 

	posAdd1(:,nu) = posAdd1(:,nu) + 1;
    posAdd2(:,nu) = posAdd2(:,nu) - 1;
        
    posP = [posP; posAdd1; posAdd2];
    
end

posP = posP(size(pos,1)+1:end,:);

end