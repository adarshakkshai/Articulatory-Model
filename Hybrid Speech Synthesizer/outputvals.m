function outputvals()
%% Inputs 
Xkk = [floor(5*randn(1,200));floor(5*randn(1,200))]';
Yk = [floor(5*randn(1,200));floor(5*randn(1,200))]';
%% Initialization matrices

Ak = [1];
Bk = [1];
Ck = [1];
Dk = [0];
Uk = [0];
%% Errors and their covariances

Qk = [1];
Rk = [1];
Pk = [1];
Pkk = [1];
Gk = [1];
I = eye(1,1);
%% Params

iterations = 100;
%%
    for i = 1:iterations
       %% Measurement
       Mk = Pkk*Ck'*inv(Rk + Ck*Pkk*Ck');
       Xk = Xkk + Mk*(Yk - Ck*Xkk - Dk*Uk);
       Pk = (I - Mk*C)*Pkk;
       %% Update
       Xkk = Ak*Xk + Bk*Uk;
       Pkk = Ak*Pk*Ak' + Gk*Qk*Gk';
       Output(:,i) = Xkk(:,1); 
    end
    %%
           plot(Xkk,'ro-'); 
%%
figure
axis tight;
axes;
vidObj = VideoWriter('trajectory.avi');
open(vidObj);
for i=1:iterations
   pause(1);
   plot(Output(:,i),'o-','markerfacecolor','r','markersize',5);
   currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
end
close(vidObj);
end