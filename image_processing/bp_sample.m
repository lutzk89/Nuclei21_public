function h = bp_sample(lf,hf,SampleSize)

[f1,f2] = freqspace(SampleSize,'meshgrid');
Hd = ones(SampleSize); 
r = sqrt(f1.^2 + f2.^2);
Hd((r<lf)|(r>hf)) = 0;

h = MYfwind1(Hd,MYhamming(SampleSize));