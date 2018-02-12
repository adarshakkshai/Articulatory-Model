%% Synthesis
ytil = liljencrantglot()
rtil = rosenbergglot()
data = load('areafnI.mat');
omega = data.sam(1:44)';
td = 500e-3;
fs = 146;
N = 44;
f1 = 2;
f2 = 3;
phi1 = 1:N;
phi2 = 1:N;
fmax = 8000;
[V,A,C] = AreaFunctionTubeTalker(omega,f1,f2,phi1,phi2,fs,td,N);
%% Obtaining input Impedances
for n = 1:td*fs
    for i = 1:fmax
        [Zin,Pout,Volv,Kn,Zl] = VocalSynthesisV1(A(:,n),i);
        Zomega(i,n) = Zin;
        Pt(:,i,n) = Pout';
        Ut(:,i,n) = Volv';
    end
end
%% Sound source convolution (glottal impulses still need to be generated)
k = [];
for i = 1:td*fs-1
    hj = 8000.*(ifft(ifftshift(fftshift(Zomega(:,i))),'symmetric'));
    k = [k conv(ytil( 1+(i-1)*55:i*55),hj(1:300),'same')];
end

%% Does de-emphasis improve hearability?
b = [1 0.95];
kb = filtfilt(b,1,k);
%%
plot((0:1/8000:0.495-1/8000),k); title('VCV Output'); xlabel('Time(s)'); ylabel('Speech Output');
soundsc(k,8192);
%%
plot((0:1/8000:0.495-1/8000),kb); title('VCV Output'); xlabel('Time(s)'); ylabel('Speech Output');
soundsc(kb,8192);