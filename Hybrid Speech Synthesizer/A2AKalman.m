function A2AKalman(x)

%%
fmax = 8000; % FFT bins required for vocal synthesis
n = 8000; % FFT Length
V = load('areafnI.mat'); % The area function
V = V.sam(1:44);
H = (conj(dftmtx(2*n))/2*n); % Measurement matrix
%% Estimate response

for i = 1:fmax
    [Zin,Pout,Volv,Kn,Zl] = VocalSynthesisV1(V(:),i);
    Zomega(i) = Zin;
    Pt(:,i) = Pout';
    Ut(:,i) = Volv';
end
%% State-Model

IT = 100;
qt = randn(1);
%Xkk = eye(n,n)*Zomega; % Assume no state noise
EE = [Zomega fliplr(Zomega)];
Yk = real(H*EE')'; % Assume no process noise
Xkk = zeros(1,16000);
Pkk = zeros(1,16000);
Xkk(1) = 10;
Pkk(1) = 1;
r = randn(1);

%%
for n = 1:IT
for k = 1:16000
    %% Update 
    K = Pkk(k) * inv(Pkk(k)+r);
    Xk = Xkk(k) + K * (Yk(k) - Xkk(k));
    Pk = (1 - K)*Pkk(k);
    %% Predict
    Xkk(k) = Xk;
    Pkk(k) = Pk + qt;
end
end

end
%% 