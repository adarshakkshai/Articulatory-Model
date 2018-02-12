function ug = glottalflow(p1,ps,png,fs,N)
% p1 is the pressure downstream of the glottis
% ps is the subglottal pressure
% png is the noise pressure generated for fricatives or aspirations
% N - Total length of signal Fourier Transform
rho = 1.147e-3; %Density of air (cm^3/s)
d1 =;
d2 =;
lg =;
mu=;

Ag01 = ; % Glottal rest area (cm^2)
Ag02 = Ag01; % Glottal rest area (cm^2), It's the same as Ag01

for n = 1:N
%% Quasi-Stationary Inductance
Ltot = rho.*(d1./Ag1(n) + d2./Ag2(n));

%% Quasi-Stationary Resistance
Rtot = 0.5.*rho( 0.37./(Ag2(n).^2) + (1 - 2*(Ag2(n)/l(n)).*(1 - Ag2(n)./l(n)))./(Ag2(n).^2)).*abs(ug(n)) + 12.*mu.*lg.^2.*(d1/Ag1(n).^3 + d2/Ag2(n).^3);

%% Glottal areas
Ag1(n) = Ag01(n) + 2*lg*x1(n);
Ag2(n) = Ag02(n) + 2*lg*x2(n);
%% Calculation of the glottal volume velocity

ug(n) = ((ps(n) - p1(n) - png(n)) + Ltot(n)*ug(n-1)*fs)./(Rtot(n) + Ltot(n)*fs);
end
end