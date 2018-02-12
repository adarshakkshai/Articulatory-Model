function A2AKalman2(x)

%%
fmax = 8000; % FFT bins required for vocal synthesis
n = 4000; % FFT Length
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
Xkk = 10*ones(n,1);
Pkk = 1*ones(n,1);
r = randn(1);
R = randn(16000,16000)
%%
for n = 1:IT
    %% Update 
    K = Pkk * H'*inv(H*Pkk*H'+ R);
    Xk = Xkk + K * (Yk - Xkk);
    Pk = (1 - K*H)*Pkk;
    %% Predict
    Xkk = Xk;
    Pkk = Pk + qt;
end

end
%% 