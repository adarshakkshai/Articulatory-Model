function rtil = rosenbergglot()
%% Different Models
f0=100;                % frequency in Hz
t=linspace(0,0.1,8000); % time axis
u=glotros(0,t*f0,ones(1,length(t)) )      % glottal waveform using default parameters
rtil = u;
plot(t,rtil);
end