clear all;

%% Area Function
t1 = 0.1:1/85:0.45;
t2 = 0.5:1/35:0.88;
omega = [sin(2*pi*t2) sin(2*pi*t1)];
data = load('areafnI.mat');
omega = data.sam(1:44)';

fs = 146; % Sampling rate (Samples)
td = 500e-3; % Duration (s)
f1 = 2; % Frequency of Mode 1
f2 = 3; % Frequency of Mode 2
phi1 = 1; % Amplitude of Mode 1
phi2 = 1; % Amplitude of Mode 2

t = 0:1/fs:td - 1/fs;
%omega = flipud(omega);
mode1 = phi1'.*cos(2*pi*f1*t);
mode2 = phi2'.*sin(2*pi*f2*t);

% V-Substrate
for n = 1:td*fs
    V(:,n) = pi/4*(omega + mode1(n) + mode2(n)).^2;        
end

% C-Substrate (No Acoustic coupling)
x =  -5:0.137:5;
y = -4:0.395685:13.4;
[X,Y] = meshgrid(x,y);
varx = 1;
vary = 1;
G = 1;
C = 1-(1/sqrt(2*pi*varx*vary))*exp(-(X.^2 + Y.^2)/(2*varx*vary));
%surf(C);

% Complete Area Function
for n = 1:td*fs
    A(:,n) = V(:,n).*C(:,n);
end

%% Synthesis of the sound (Remember always pass row)
St = [];
for n = 1:td*fs
    for i = 1:8000
     [Zin,Pout,Volv,Kn,Zl] = VocalSynthesisV1(V(:,n),i);
     Zintotal(i,n) = Zin;
    end
    % St = [St Pout(:,n)'];
    %end
end
figure;
plot(abs(Zintotal));
%figure;
%plot(flipud(omega)')
%%
hj = Zintotal;
% max(abs(hj))
% hj(:,8:12) = 0;
% hj(:,63:67) = 0;
% hj(:,58) = 0;
% hj(:,17) = 0;
% hj(:,51) = 0;
% hj(:,24) = 0;
hj(:,50:52) = 0;
hj(:,23:25) = 0;
max(abs(hj))

yo = [];
for i = 1:73
    ht = (ifft(ifftshift(fftshift(hj(:,i))),'symmetric'));
    ou = conv(ug,ht(1:512),'same');
    yo = [yo ou];
end
plot(yo);
soundsc(yo,8000);
plot(abs(hj));

%%
figure;
h = surf(t*1000,(1:44)*0.396825,V);
set(h,'LineStyle','none');
%%
figure;
g = surf(t*1000,(1:44)*0.396825,A); axis xy;
set(g,'LineStyle','none');
close all;
%% Kalman Using Le-Talker for observation
%A = []; % These values go with the system state, we fill in values of the parameters distributed in f0
% Assume 
%B = []; % These values go with the control state, the control parameters are closely linked to the supra glottal and sub glottal thresholds 
%C = []; % Measurement indices with respect to requirement
%wN = 0; % Noise of source
%wT = 0; % Noise of measurement

