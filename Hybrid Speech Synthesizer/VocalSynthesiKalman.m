function Pout = Vocalsynthesis(arfn)
%% Constants
length = 0.396825; % Tube length (cm)
tau = 35000.0;% Speed of sound (m/s)
rho = 1.14e-3; % Density of air (kg/m^3)
b = 3;       % Radius of the lip (cm)
fw = 15;     % Mechanical resonance frequency of wall (Hz)
r = 408;     % Wall resistance to mass (rad/s)
ft = 200;    % Lowest resonance frequency of vocal tract (Both ends closed) (Hz)
q = 4;       % Correction for thermal conductivity and viscosity (rad/s)
f = 900;     % Frequency (Hz)
Pin = 200;   % Initial glottal pressure (kg/ms^2)
Vin = 7;    % Initital glottal velocity (m/s)
n = 44;      % Total number of sections (Unitless)
%% Transmission-Line Model for Pressure output

% Impedance
Am = pi*b^2;
Zm = rho*tau/Am;
R = 128*Zm/(9*pi^2);
L = 8*b*Zm/(3*pi*tau);
Zl = (R*L^2*(2*pi*f)^2 + 1i*2*pi*f*R^2*L)/(R^2 + (2*pi*f)^2*L^2);

% Wave guide parameters
al = sqrt(1i*2*pi*f*q);
bt = ( (1i*2*pi*f*(2*pi*ft)^2)/((1i*2*pi*f + r)*(1i*2*pi*f) + (2*pi*fw)^2) )+al;
gam = sqrt( (r + 1i*2*pi*f)/(bt + 1i*2*pi*f));
sig = gam*(bt + 1i*2*pi*f);
A = cosh(sig*length/tau);
D = cosh(sig*length/tau);

% Initial Pressures subglottal
Pout = [Pin];
Volv = [Vin];

% Output Pressures
for l = 2:n
        B = -rho*tau*gam*sinh(sig*length/tau)/arfn(l-1);
        C = -arfn(l-1)*sinh(sig*length/tau)/(rho*tau*gam); 
        Pout(l) = A*Pout(l-1) + B*Volv(l-1);
        Volv(l) = C*Pout(l-1) + D*Volv(l-1);
end


% Lip Impedance
Pout = Pout*Zl;

end