function w = hammings(M)
w = .54 - .46*cos(2*pi*(0:M-1)'/(M-1));
