function reihentest(y,x)
set(0,'RecursionLimit',1000)
global imbin image siz num sum x_pivot y_pivot
x_start=x;
x_end=x;

% search left:
while x_start>1 & imbin(y,x_start-1)==1
    x_start=x_start-1;
    pixelreport(y,x_start)
end

% search right:
while x_end<siz(2) & imbin(y,x_end+1)==1
    x_end=x_end+1;
    pixelreport(y,x_end)
end

% search down:
if y>1
    for xs=x_start:x_end
        if imbin(y-1,xs)==1
            pixelreport(y-1,xs)
            reihentest(y-1,xs)
        end
    end
end

% search up:
if y<siz(1)
    for xs=x_start:x_end
        if imbin(y+1,xs)==1
            pixelreport(y+1,xs)
            reihentest(y+1,xs)
        end
    end
end