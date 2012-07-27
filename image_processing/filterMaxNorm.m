Pow    = 100.;
preFac = 100.;
invPreFac = 1./preFac;
invPow    = 1./Pow;

im100pre = max(im - thresh(im, thresholdSkew), 0);

if FigNums(6)>0
    figure(FigNums(6));subplot(2,2,1);imagesc(im100pre);colorbar
end
im100    = (preFac*im100pre).^Pow;
im100filter = (filter2(circle(circleSize100filter),im100).^invPow).*invPreFac;
% 
if FigNums(6)>0
    figure(FigNums(6));subplot(2,2,2);imagesc(im100filter);colorbar
end

im100 = im100pre-im100filter*subtract100fac;
im100(im100<0) = 0;

if FigNums(6)>0
    figure(FigNums(6));subplot(2,2,3);imagesc(im100);colorbar
end

%figure(53);imagesc(im100>0.01)
