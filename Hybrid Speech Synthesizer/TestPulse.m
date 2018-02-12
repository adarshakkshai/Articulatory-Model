clear all;
%% Glottal Source

t = 0:1/44100:1- 1/44100;
h = fliplr(0:1/22050:1);
u = zeros(1,44101);
u(22050:end-1) = h.*cos(2*pi*t(1:22051));%+0.2925;
%plot(u)
u(1:22050) = cos(2*pi*t(22050:end-1));
%plot(u(11027:33075));
u = u(11027:33075);
u = [ u zeros(1,22051)];
%plot(u)
ug = u;

%%
t = 0 : 1/8000 : 1 - 1/8000;
d = 5e-3 : 1/50  :  1 - 1/50;
y = pulstran(t,d,'rectpuls',10e-3); 
plot(t,y); hold on;
ug = y;
 FS = 8000;
 F0 = 100;
 t = linspace(0,1,FS);
 y = zeros(size(t));
 y(1:FS/F0:end)=1;
 y(2:FS/F0:end)=1;
 y(3:FS/F0:end)=1;
ug = y;
plot(t,ug);
hk = (0.95).^((1:128));
hk = (conv(hk,fliplr(hk)))';
%% Area Function
t1 = 0.1:1/85:0.45;
t2 = 0.5:1/35:0.88;
omega = [sin(2*pi*t2) sin(2*pi*t1)];
data = load('areafnI.mat');
omega = data.sam(1:44)'; 

%% Input characteristic impedance of first section
rho = 1.147e-3;  % Density of air (g/cm^3)
tau = 35000.00;    % Speed of sound (m/s)
length = 0.396825;  % Length of each section (cm)
Zo = rho*tau/length;

%% Obtaining input Impedances
for i = 1:8000
    [Zin,Pout,Volv,Kn,Zl] = VocalSynthesisV1(omega,i);
    Zomega(i) = Zin;
    Pouttotal(:,i) = Pout';
end
%% Plots for report
plot(20*log(abs(Zomega)),'b');
axis([0 3000 0 180])
set(gca,...
'Units','normalized',...
'YTick',0:10:180,...
'XTick',0:200:length(Zomega),...
'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times');

ylabel({'Log Magnitude of Impedance'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',9,...
'FontName','Times');
xlabel('Frequency(Hz)',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',7,...
'FontName','Times');

legend({'$\log(Z_{in})$'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',7,...
'FontName','Times',...
'Location','NorthEast');
title('Vocal Tract Input Impedance',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',7,...
'FontName','Times')
print -depsc2 ZinPlot.eps
%%
clear hj;
hj = (ifft(ifftshift(fftshift(Zomega)),'symmetric'));
figure;
plot((hj));
gh(1:25) = 0;
gh(26:50) = 0;
gh(51:75) = 1;
gh(76:100) = 1 ;
figure;
plot(conv(gh,hj(1:256)));
plot([0:1/8000:0.05-1/8000],conv(ug(1:400),hj(1:100),'same'));
soundsc(conv(ug,hj,'same'),8000);
%% Adding Lip imepdance to Tract bank
for i = 0:43
    Zin(i+1) = (Kn(2,i*2+2)*Zl - Kn(1,i*2+2))/(Kn(1,i*2+1) - Kn(2,i*2+1));
end
%% Assuming delay banks

for i = 0:43
    b(i+1) = (Kn(2,i*2+2)*Zl - Kn(1,i*2+2));
end

for i = 0:43
    a(i+1) = (Kn(1,i*2+1) - Kn(2,i*2+1));
end


%% Input Reflectance
N = 44100;
for i = 1:44
Rin(i) = (Zin(i) - Zo)/(Zin(i) + Zo); % Frequency domain
end
rin = ifft(Rin,N,'symmetric');  % Time domain, with N Samples

%%
rin = rin(1:8001);
zt = ifft(Zin,250,'symmetric');
%% The convolution to obtain supraglottal pressure
pl = zeros(1,N); % Initital supraglottal pressure

for n = 1:N
   first = ((1+rin(1))/(1-rin(1)))*Zo*ug(n);
   rintemp = rin(2:end);
   second = conv(rintemp,pl);
   third = Zo*conv(rintemp,ug);
   fourth = second(n) + third(n);
   fifth = (1/(1 - rin(1)))*fourth;
   pl(n) = first + fifth;
   n
end
%% Alternatively
plout = [pl(1:11126)];
plsound = [];
for i = 1:200
plsound = [plsound plout zeros(1,11025)];
end
plout = conv(zt,ug);