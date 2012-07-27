function moments=findFrq(im, bins, plotfac)
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
rindices = 1+floor(bins*Mfrq);
%figure(21);imagesc(rindices);colorbar
imFFTmean = accumarray(rindices(:),imFFT(:))./accumarray(rindices(:),1);
frIndices=(0:length(imFFTmean)-1)/bins;
figure(13);plot(frIndices(3:round(bins*plotfac)),(imFFTmean(3:round(bins*plotfac))));
normal=sum(sum(imFFT));
moments=sqrt(sum(sum(imFFT.*Mfrq2))/normal);
