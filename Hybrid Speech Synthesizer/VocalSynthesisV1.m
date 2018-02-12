function [Zin,Pout,Volv,Kn,Zl] = VocalSynthesisV1(arfn,f)
%% Constants (Not sure about the values of Pin and Vin)
length = 0.396825; % Tube length (cm)
tau = 35000.0;% Speed of sound (cm/s)
rho = 1.14e-3; % Density of air (g/cm^3)
%b = 0.3554;       % Radius of the lip (cm)
Am = arfn(end);    % Area near termination of the lips
fw = 15;     % Mechanical resonance frequency of wall (Hz)
r = 408;     % Wall resistance to mass (rad/s)
ft = 200;    % Lowest resonance frequency of vocal tract (Both ends closed) (Hz)
q = 4;       % Correction for thermal conductivity and viscosity (rad/s)
Pin = 0.0007;   % Initial glottal pressure (dyn/cm^2)
Vin = 1;    % Initital glottal volume velocity (cm^3/s)
n = 44;      % Total number of sections (Unitless)
%% Transmission-Line Model for Pressure output

%% Radiation Impedance
%%
%Am = pi*b^2;
%Zm = rho*tau/Am; % Characteristic acoustic impedance at the end of lips
%%
b = sqrt(Am/pi);
Zm = rho*tau/Am; % Characteristic acoustic impedance at the end of lips
R = 128*Zm/(9*pi^2);
L = 8*b*Zm/(3*pi*tau);
Zl = (R*L^2*(2*pi*f)^2 + 1i*2*pi*f*R^2*L)/(R^2 + (2*pi*f)^2*L^2);

%% Wave guide parameters
al = sqrt(1i*2*pi*f*q);
bt = ( (1i*2*pi*f*(2*pi*ft)^2)/((1i*2*pi*f + r)*(1i*2*pi*f) + (2*pi*fw)^2) )+al;
gam = sqrt( (r + 1i*2*pi*f)/(bt + 1i*2*pi*f));
sig = gam*(bt + 1i*2*pi*f);
A = cosh(sig*length/tau);
D = cosh(sig*length/tau);

%% Initial Pressures subglottal
Pout = [Pin];
Volv = [Vin];

%% Frequency Dependent Matrix 
% This Calculation of B and C is  from Story (2017) TubeTalker Model
% Many ways to compute these two elements B and C

B = -rho*tau*gam*sinh(sig*length/tau)/arfn(1);
C = -arfn(1)*sinh(sig*length/tau)/(rho*tau*gam);
Pout(1) = A*Pin + B*Vin;
Volv(1) = C*Pin + D*Vin;
Kn = [A B;
      C D];
%% Output Pressures and velocities
for l = 2:n
        B = -rho*tau*gam*sinh(sig*length/tau)/arfn(l);
        C = -arfn(l)*sinh(sig*length/tau)/(rho*tau*gam); 
        Pout(l) = A*Pout(l-1) + B*Volv(l-1);
        Volv(l) = C*Pout(l-1) + D*Volv(l-1);
        K = [A B;
             C D];
        Kn = [K*Kn];
end

%% Vocal tract output pressure at lips
%Pout = Volv(end)*Zl;

%% Lip Impedance
Zin = (Kn(2,end)*Zl - Kn(1,end))/(Kn(1,end-1) - Kn(2,end-1)*Zl);

end