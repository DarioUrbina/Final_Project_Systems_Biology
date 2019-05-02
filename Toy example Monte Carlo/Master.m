%% Differential equation details
tstep = 0.01;
stopTime= 2;
resolution = (stopTime / tstep)+1;
tspan = [0:tstep:stopTime]; 
n_diffeqn = 1;
options = odeset('RelTol',1e-6); 
initvalue = zeros(n_diffeqn,1);
initvalue(1,1) = 6; %Y
iterations = 100;
%% Parameters
% params(1,1) = 1; 
% params(2,1) = 11;
% params(3,1) = 8; 
% params(4,1) = 0.2; 
%% First parameter
params(1,1:resolution) = 1; 
%% Third parameter
params(3,1:resolution) = 8; 
%% Fourth parameter
params(4,1:resolution) = 0.2; 
%% Figure 1
figure
maxY = 10; 
axis([0 stopTime 0  maxY]);
%% Results vector
results = zeros(resolution,iterations);

%%
for iteration=1:iterations 
%% Second parameter
a_1 = 7 ;
a_2 = 13 ;
params_2 = a_1 + (a_2-a_1).*rand(1,resolution);
params(2,1:resolution) = params_2;
%%
    %for c=1:size(tspan,2)
    c=1;
    [tsim,results(:,iteration)] = ode15s(@Function_B,tspan,initvalue,options,params(:,c));
    %end
    %IMPORTANT: A random value for a parameter is used for the whole time span of the
    %experiment. The parameter value does't change through the time span of
    %the experiment. Initially I was having that value changig for each
    %step intime, but I was taking only the last step! 
    
plot(tsim,results(:,iteration),'LineWidth',.7);
hold on

end

%legend ('A','B')
xlabel('time(s)')
ylabel('concentration (M)')
title ('Transient part')
results;

%% Figure 2
top = max(results(3,:))
top_idx = find(results(3,:) == top)
bottom = min(results(3,:))
bottom_idx = find(results(3,:) == bottom)

figure

plot(tsim,results(:,top_idx),'LineWidth',.7);
hold on
plot(tsim,results(:,bottom_idx),'LineWidth',.7);

% hold on

% Plot results
% figure(3);
% figure
% maxY = 10; 
% plot(tsim,results,'LineWidth',4);
% axis([0 stopTime 0  maxY]);
% legend ('A','B')
% xlabel('time(s)')
% ylabel('concentration (M)')
% title ('Transient part')
% Plot results
% figure(4);
% plot(tsim4,results4,'LineWidth',4);
% axis([0 stopTime_2 0  maxY]);
% legend ('A','B')
% xlabel('time(s)')
% ylabel('concentration (M)')
% title ('Steady-State part')

% %%
% params(1,1) = 1; 
% params(2,1) = 11;
% params(3,1) = 8; 
% params(4,1) = 0.2; 
% % Set details
% tstep = 0.2;
% stopTime= 2;
% tspan = [0:tstep:stopTime]; 
% 
% %for c = 1:s
% %size(tspan)
% 
% % for c=1:size(tspan,2)
% %     c
% % end
% 
% maxY = 10; 
% n_diffeqn = 1;
% options = odeset('RelTol',1e-6); %options for ODE
% 
% % Specify initial values
% initvalue = zeros(n_diffeqn,1);
% initvalue(1,1) = 6; %Y
% 
% % Refer to file with ODEs and call solver
% [tsim,results] = ode15s(@Function_B,tspan,initvalue,options,params)
% % Refer to file with ODEs and call solver
% 
% % Plot results
% figure(3);
% plot(tsim,results,'LineWidth',4);
% axis([0 stopTime 0  maxY]);
% legend ('A','B')
% xlabel('time(s)')
% ylabel('concentration (M)')
% title ('Transient part')
% % Plot results
% % figure(4);
% % plot(tsim4,results4,'LineWidth',4);
% % axis([0 stopTime_2 0  maxY]);
% % legend ('A','B')
% % xlabel('time(s)')
% % ylabel('concentration (M)')
% % title ('Steady-State part')
% 
%  r = a + (b-a).*rand(N,1)
%  a=20
%  b=150