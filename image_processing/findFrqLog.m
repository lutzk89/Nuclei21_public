function moments=findFrqLog(im, bins, skew, LocalMaxKernel, plotfac, FigNum)
if nargin<6
    FigNum = -1;
end
if nargin<4
    LocalMaxKernel = -1;
end
    
siz=size(im);
%siz2=round(siz/2.);
imFFT=abs(fftshift(fft2(im)));%.^2;
%imFFTpol=cart2pol(imFFT);
%figure(100);image(imFFT*.01);
[frq1 frq2] = freqspace(siz,'meshgrid');

% frq1=0:siz(1)-1;
% frq2=0:siz(2)-1;
% frq1(siz2(1)+1:siz(1))=frq1(siz2(1)+1:siz(1))-siz(1);
% frq2(siz2(2)+1:siz(2))=frq2(siz2(2)+1:siz(2))-siz(2);
% %cart2pol(frq1,frq2)
% %sum(frq1)
% frq1=frq1';
Mfrq2 = frq1.^2 + frq2.^2;
Mfrq  = sqrt(Mfrq2);
MfrqFunc   = max(log10(Mfrq)+3,0);

%Maxfrq     = max(max(Mfrq));
MaxfrqFunc = log10(1)+3;

rindices = 1+floor(bins*MfrqFunc/MaxfrqFunc);
r0index  = 1+floor(bins*(log10(1/120)+3)/MaxfrqFunc);
r1index  = 1+bins;
%rindices = 1+floor(bins*Mfrq);
%figure(21);imagesc(rindices);colorbar
imFFTmean = accumarray(rindices(:),imFFT(:))./accumarray(rindices(:),1);

if skew~=0
    tiltLine = 1 + skew*((1-r0index):(length(imFFTmean)-r0index))/length(imFFTmean);
    imFFTmean = imFFTmean.*tiltLine';
end
%frIndices=(0:length(imFFTmean)-1)/bins;
%figure(13);plot(frIndices(1:round(bins*plotfac)),(imFFTmean(1:round(bins*plotfac))));
imFFTfiltered = conv(imFFTmean,[.25 .5 .25],'same');
imFFTfilteredRange = imFFTfiltered(r0index:r1index);

if length(LocalMaxKernel)==1 && LocalMaxKernel<=0
    [m, i] = max(imFFTfilteredRange);
else
    [m, i] = localMax(imFFTfilteredRange, LocalMaxKernel);
end
%[m, i] = max(imFFTmean(r0index:r1index));

i = i+r0index-1;
if FigNum>0
    figure(FigNum);
    plot(imFFTmean(10:bins*plotfac),'g');
    hold on;
    plot(imFFTfiltered(10:bins*plotfac))
    plot(i-9,m,'ro');
    hold off;
end
moments = 2/ifunc(i,bins);
%normal=sum(sum(imFFT));
%moments=sqrt(sum(sum(imFFT.*Mfrq2))/normal);
