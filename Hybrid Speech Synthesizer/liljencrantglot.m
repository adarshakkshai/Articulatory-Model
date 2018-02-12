function ytil = liljencrantglot()
%% Different Models
f0=100;                % frequency in Hz
t=linspace(0,0.5,4000); % time axis
u=glotlf(0,t*f0)       % glottal waveform using default parameters
y = fir1(100,0.5,'low');
ytil = filtfilt(y,1,u);
plot(t,ytil);
end
