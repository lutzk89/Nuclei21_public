function pixelreport(y,x)
global imbin image siz num sum x_pivot y_pivot

imbin(y,x)=2;
x_pivot = x_pivot + x*double(image(y,x));
y_pivot = y_pivot + y*double(image(y,x));

num=num+1;
sum=sum+double(image(y,x));
