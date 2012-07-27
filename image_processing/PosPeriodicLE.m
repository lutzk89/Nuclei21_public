function posP = PosPeriodicLE( pos, LEshift, rim1 )

dim = size(LEshift,1);
L   = zeros(dim,1);

for nu=1:dim
    L(nu) = LEshift(nu,nu);
end

disp(L);

rim2 = L-rim1;

posP = pos;
for nu=1:dim %cycle through all boundaries
    posAdd1 = posP(posP(:,nu)<rim1(nu),:); 
    posAdd2 = posP(posP(:,nu)>rim2(nu),:); 
    for mu=1:dim %cycle through all directions
        posAdd1(:,mu) = posAdd1(:,mu) + LEshift(nu,mu);
        posAdd2(:,mu) = posAdd2(:,mu) - LEshift(nu,mu);
        
        %if mu>nu
        %    posAdd1(:,mu) = mod( posAdd1(:,mu) , L(:,mu) );
        %    posAdd2(:,mu) = mod( posAdd2(:,mu) , L(:,mu) );
        %end
    end
    
        
    posP = [posP; posAdd1; posAdd2];
    
end

end