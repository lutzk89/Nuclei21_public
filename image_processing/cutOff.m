function out=cutOff(x,cutoff0,cutoff1)
out=x;
out(x<cutoff0)=cutoff0;
out(x>cutoff1)=cutoff1;