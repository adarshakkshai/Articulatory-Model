function [E] = EventFunction(N1,Npk,N2,td,fs)
%N1 - Start sample of deflection
%Npk - Peak sample of deflection
%N2 - Stop sample of deflection
%td - Time duration of the utterance
%fs - Sampling rate of the utterance
%As per Story et.al in their 2017 paper

for n = 1:td*fs
    if ( n <=N1)
        E(n) = 0;
    elseif (n>N1 && n<= Npk)
        E(n) = 0.5.*(1 - cos(pi.*(n-N1)./(Npk-N1)));
    elseif (n>Npk && n <= N2)
        E(n) = 0.5.*(1 + cos(pi.*(n-Npk)./(N2-Npk)));
    elseif (n>N2 && n<=td*fs)
        E(n) = 0;
    end
end
end