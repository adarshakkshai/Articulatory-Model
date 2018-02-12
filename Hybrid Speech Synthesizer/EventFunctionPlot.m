%% Plotting the event function
E = EventFunction(5,30,55,500e-3,146);
plot((0:1/146:500e-3-1/146).*1e3,E); hold on;
k = nan*zeros(1,73);
k(1,5) = 1;
k(1,55) = 1;
stem((0:1/146:500e-3-1/146).*1e3,k,'x');
%%