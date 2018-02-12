function [V,A,C] = AreaFunctionTubeTalker(Iarea,f1,f2,phi1,phi2,fsa,td,N)
%Iarea - Initial Area function to be introduced into the TubeTalker Model (cm^2)
%f1 - Modulation frequency of Q1 (Hz)
%f2 - Modulation frequency of Q2 (Hz)
%phi1 - Amplitude ramps of Q1 (Const./Vector)
%phi2 - Amplitude ramps of Q2 (Const./Vector)
%fs - sampling frequency of the area function slices (Hz)
%td - time duration of the area function (s)
%N - Number of sections of the vocal tract(Typically 44)
%vtlength - Total length of the assumed VT (Typically 17.4603cm)
if size(Iarea,1) == 1
    Iarea = Iarea.';
end

if size(phi1,1) == 1
    phi1 = phi1.';
end

if size(phi2,1) == 1
    phi2 = phi2.';
end

if (size(phi1,1) ~= N) && (size(phi2,1) ~= N)
    error('The amplitude vector needs to be equal to the length of sections');
end
vtlength = 17.4603; % VT length (cm)
dl = round(vtlength/N,6);
fmax = 8000;
%% Calculate the Cosine modes

t = 0 : 1/fsa : td - 1/fsa;
mode1 = phi1.*cos(2*pi*f1*t); % Assume In-phase 
mode2 = phi2.*sin(2*pi*f2*t); % Assume pi/2 phase shift
%mode1 = feval('simplesine',f1,phi1,t);
%mode2 = feval('simplecosine',f2,phi2,t);
%% Vowel Substrate 

for n = 1:td*fsa
     V(:,n) = Iarea;
    %Using this would require the Iarea to be a diameter function and not so
    %much of a area function with a radius(Beware of this error)
    %V(:,n) = pi/4*(Iarea).^2 + mode1(:,n) + mode2(:,n)).^2;        
end
%% Consonant Substrate (No Acoustic coupling)
%Random Consonant Substrate
%x =  -5:0.137:5;
%y = -4:0.395685:13.4;
%[X,Y] = meshgrid(x,y);
%varx = 1;
%vary = 1;
%C = 1-(1/sqrt(2*pi*varx*vary))*exp(-(X.^2 + Y.^2)/(2*varx*vary));
%% %% Consonant Substrate (Practical Method)
%Acoustically coupled to the area function
repi = 6; %Numerator of substrate
iepi = 5; %Denominator of substrate
N1 = 5; %Start peak
Npk = 30; %Primary peak
N2 = 55; %(Secondary die down peak)
mu = 1; % (Random mu parameter as described in Story et.al. 2017)
rlimit = 3; % Resonance limit ( only first 3 resonances are considered)
delta = [-1 -1 0]; %deflections in the consonants
[C] = Csubstrate(td,fsa,V,dl,fmax,N,repi,iepi,delta,N1,Npk,N2,mu,rlimit);
%% Complete Area Function
% Why is the scaling so low? It's supposed to be above 0 and between [0,2]
%It's actually between [-1 1], buggy code? Recheck again
C = 1+C;
for n = 1:td*fsa
    A(:,n) = V(:,n).*C(:,n);
end
end

