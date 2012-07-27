function Hist=dichtehist(datei,bereich,L,N,i); 
% L:laenge der box; 
% N: Anzahl Bins in jeder Dimension;

%n = dlmread(datei, '', bereich);

Nteilchen = 1000000;
% zurückholen in die box
xd = rand(Nteilchen ,1)*L;%mod(n(:,1),L); %mod L
yd = rand(Nteilchen ,1)*L;%mod(n(:,2),L);
zd = rand(Nteilchen ,1)*L;%mod(n(:,3),L);

% % Linie zerlegen bzw. Volumen zerlegen
% xi = linspace(0,L,N);
% yi = linspace(0,L,N);
% zi = linspace(0,L,N);
% 
% % Interpolieren
% xr = interp1(xi,0.5:numel(xi)-0.5,xd,'nearest');
% yr = interp1(yi,0.5:numel(yi)-0.5,yd,'nearest');
% zr = interp1(zi,0.5:numel(zi)-0.5,zd,'nearest');

xr = floor(xd/L*N)+1;
yr = floor(yd/L*N)+1;
zr = floor(zd/L*N)+1;


% Zusammenzählen
z = accumarray([xr yr zr], 1, [N N N]);

% Mache aus 3d-array z einen Linienvektor
zn = reshape (z,numel(z),1);


% Histogrammzählen
histline = 0:1:max(zn);
Hist = histc(zn,histline);


figure(i);

plot(histline,Hist,'-');

hold on
mu=Nteilchen/N^3;
HistPois = N^3* exp(-mu)*(mu.^histline)./factorial(histline);
plot(histline,HistPois,'r-');
hold off
%cd /usr/scratch2/fiege/