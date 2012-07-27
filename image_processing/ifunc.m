function out = ifunc( i, bins )
out = 10.^((i-1) / bins * (log10(1)+3) - 3);
end

