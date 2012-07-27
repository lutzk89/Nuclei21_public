function h = fsamp2(hd) 
%FSAMP2 Design 2-D FIR filter using frequency sampling. 
%   FSAMP2 designs two-dimensional FIR filters based on a desired 
%   two-dimensional frequency response sampled at points on the 
%   Cartesian plane. 
% 
%   H = FSAMP2(HD) designs a two-dimensional FIR filter with 
%   frequency response HD, and returns the filter coefficients in 
%   matrix H. (FSAMP2 returns H as a computational molecule, 
%   which is the appropriate form to use with FILTER2.) The 
%   filter H has a frequency response that passes through points 
%   in HD. If HD is M-by-N then H is also M-by-N. 
% 
%   HD is a matrix containing the desired frequency response 
%   sampled at equally spaced points between -1.0 and 1.0 along 
%   the x and y frequency axes, where 1.0 corresponds to half the 
%   sampling frequency, or pi radians. For accurate results, use 
%   frequency points returned by FREQSPACE to create HD. 
% 
%   H = FSAMP2(F1,F2,HD,[M N]) produces an M-by-N FIR filter by 
%   matching the filter response HD at the points in the vectors 
%   F1 and F2. The resulting filter fits the desired response as 
%   closely as possible in the least squares sense. For best 
%   results, there must be at least M*N desired frequency 
%   points. FSAMP2 issues a warning if you specify fewer than M*N 
%   points. 
% 
%   Class Support 
%   ------------- 
%   The input matrix HD must be numeric and nonsparse. All other  
%   inputs to FSAMP2 must be of class double. All outputs are of  
%   class double. 
% 
%   Example 
%   ------- 
%   Use FSAMP2 to design an approximately symmetric 
%   two-dimensional bandpass filter with passband between 0.1 and 
%   0.5 (normalized frequency). 
%    
%       [f1,f2] = freqspace(21,'meshgrid'); 
%       Hd = ones(size(f1)); 
%       r = sqrt(f1.^2 + f2.^2); 
%       Hd((r<0.1) | (r>0.5)) = 0; 
%       h = fsamp2(Hd); 
%       freqz2(h) 
% 
%   See also CONV2, FILTER2, FREQSPACE, FTRANS2, FWIND1, FWIND2. 
 
%   Copyright 1993-2002 The MathWorks, Inc.   
%   $Revision: 5.20 $  $Date: 2002/03/28 15:39:46 $ 
 
% Reference: Jae S. Lim, "Two Dimensional Signal and Image Processing", 
%            Prentice Hall, 1990, pages 213-217. 
 
%checknargin(1,4,nargin,mfilename); 
 
%if nargin==1, % Uniform spacing case (fast) 
%  hd = f1; 
 
hd = rot90(fftshift(rot90(hd,2)),2); % Inverse fftshift 
h = fftshift(ifft2(hd)); 
