function outAll = num2strn(numbers,n)
if nargin == 1
    n = 1;
end

outAll = [];
for j = 1:length(numbers)
    out = num2str(numbers(j));
    while length(out)<n
        out = ['0' out];
    end
    outAll = [outAll; out];
end
