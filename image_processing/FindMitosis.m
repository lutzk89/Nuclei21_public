t1=filenums(1)+5;
t2=filenums(end)-5;
array = log10(meanEdgLengths(t1:t2));
chi2  = @(Args) sum( (array(1:round(Args(2))-1) - Args(1)).^2 ) +  sum( (array(round(Args(2)):end) - Args(3)).^2 );
ArgsOut = fminsearch(chi2,[array(1),round(length(array)/2),array(end)]);

figure(15);
plot(array)
hold on;
fitArray = ArgsOut(1)*ones(length(array),1);
fitArray(round(ArgsOut(2)):end) = ArgsOut(3);
plot(fitArray,'r');
hold off;

MitosisTime = t1-1 + ArgsOut(2);
