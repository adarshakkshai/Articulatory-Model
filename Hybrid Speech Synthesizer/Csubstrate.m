function [C] = Csubstrate(td,fs,V,dl,fmax,N,repi,iepi,delta,N1,Npk,N2,mu,rlimit)
%l - sectional length (cm)
%V - Vowel area function (cm^2)
%fmax - Maximum frequency of resonance (Hz)
%td - Time duration(s)
%fs - Sampling rate (Hz)
%N - Sectional count
%rlimit - resonance limit; j =1,2,3 (first 3 resonances)
%repi - Numerator of event function
%iepi - Denomincator of event function
% Keep in mind always that the functions being generated is for the forward
% model, the reverse model will not utilize these and might need to
% prediction done in another mechanism, maybe phi's can be computed by
% direct inversion nothing more.
%% Calculate Resonances first to determine where to modulate
Pt = [];
Ut = [];
for n = 1:td*fs
    for i = 1:fmax
        [Zin,Pout,Volv,Kn,Zl] = VocalSynthesisV1(V(:,n),i);
        Zomega(i,n) = Zin;
        Pt(:,i,n) = Pout';
        Ut(:,i,n) = Volv';
    end
    [peak,location] = findpeaks(abs(Zomega(:,n)));
    %% This is not possible as the area function touches zero and causes no resonances in the given
    % section of the file
    if size(peak,1) ~= 8
        padval = 8 - size(peak,1);
        peak = padarray(peak,[padval,0],'post','symmetric');
        location = padarray(location,[padval,0],'post','symmetric');
    end
    pk{n} = peak';
    loc{n} = location';
end
%% Calculate the energy of the resonances
tau = 35000; %Speed of sound (cm^2/s)
rho = 1.147e-3; % Density of air (g/cm^3)
rlimit = 3; % Total of 3 resonances
%% Sensitivity
% Does not take differential lengths, has to be static length
for n = 1:td*fs
    for l = 1:rlimit
        Uj = Ut(:,loc{n}(l),n);
        Pj = Pt(:,loc{n}(l),n);
        KE(:,l,n) = 0.5.*rho.*dl.*(1./V(:,n)).*(Uj.*conj(Uj));
        PE(:,l,n) = 0.5.*dl.*(V(:,n)./(rho*tau.^2)).*(Pj.*conj(Pj));
        TE(:,l,n) = sum(KE(:,l,n)+PE(:,l,n));
        Sj(:,l,n) = (KE(:,l,n) - PE(:,l,n))/TE(:,l,n);
    end
end
%% Deflections
for n = 1:td*fs
    Sns = [];
    for l = 1:rlimit
        Sns = vertcat(Sns,Sj(:,l,n)');
    end
    y0(:,n) = (delta*Sns)';
end
%% Glottal effects
g = GlottalConstraint(N,repi,iepi);
% Deflections with glottal effects
%gt = repmat(g,1,N);
y = g.*y0;
%% Eventfunction 
E = EventFunction(N1,Npk,N2,td,fs);
%% C-Substrate
for n = 1:td*fs
    C(:,n) = -(mu.*E(n).*y(:,n))./min(y(:,n));
end
end