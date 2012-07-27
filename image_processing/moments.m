function moments=findFrq(im)
siz=size(im);
siz2=round(siz/2.);
imFFT=fft2(im);
frq1=1:siz(1);
frq1(siz2(1):siz(1))=frq1(siz2(1):siz(1))-siz(1);
frq2(siz2(2):siz(2))=frq2(siz2(2):siz(2))-siz(2);



moments=1
