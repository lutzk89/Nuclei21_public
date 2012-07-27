function out=num2str4(i)
pre = '';

    if i<1000 pre = [pre '0']; end
    if i<100  pre = [pre '0']; end
    if i<10   pre = [pre '0']; end
    
out = [pre num2str(i)];