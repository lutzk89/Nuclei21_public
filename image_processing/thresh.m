function threshold = thresh( I, tilt )
% Find threshold value using Otsu's method

I = I(:);
Imax = max(I);
Imin = min(I);

IbinMax = 255;
Ibin = uint8(floor( IbinMax* (I-Imin)/(Imax-Imin) ))';
%sigma2 = var(I);

range  = 0:IbinMax;
Ihist  = histc(Ibin,range);
Ihisti = Ihist.*range;

csum1i = cumsum(Ihisti);
csum1  = cumsum(Ihist);

csum2i = csum1i(IbinMax) - csum1i;
csum2  = csum1(IbinMax)  - csum1;

mu1 = csum1i./csum1;
mu2 = csum2i./csum2;

tilts = tilt*linspace(-1,1,IbinMax+1) + 1;

[~, threshold] = max(tilts.*csum1.*csum2.*(mu2 - mu1).^2);

threshold = Imin + double(threshold)/IbinMax*(Imax-Imin);
end
