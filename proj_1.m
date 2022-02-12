 clear,
close all,

lambda = 0.001; %decays/s
M = 5000; % total number of time steps
t = 0:(M-1); % time series
N0 = 3e3;
NA = zeros(M,1); 
NB = zeros(M,1); 
NA(1) = N0;

Time1=tic;
for j=2:M
    NA(j) = NA(1)*exp(-lambda*j); % number of type 1 atoms that are left
    NB(j) = NA(1) - NA(j); % number of type 2 atoms   
end
Time2=toc(Time1);
f1=figure;
figure(f1);
plot(t, [NA, NB, NA+NB]), box off
xlabel('Time (s)');
ylabel('Number of atoms');
title('Radioactive Decay (exponential)','FontSize',14);

legend({'# atoms type A', '# atoms type B','# total atoms'})
legend boxoff

NC = zeros(M,1); 
ND = zeros(M,1); 
NC(1) = N0;
Time3=tic;
for j=2:M
    n_NCtoND = nnz(rand(NC(j-1),1) <= lambda);
    NC(j) = NC(j-1) - n_NCtoND; % number of type 1 atoms that are left
    ND(j) = ND(j-1) + n_NCtoND; % number of type 2 atoms   
end
Time4=toc(Time3);
f2 = figure;
figure(f2);
plot(t, [NC, ND, NC+ND]), box off
xlabel('Time (s)');
ylabel('Number of atoms');
title('Radioactive Decay (random)','FontSize',14);
legend({'# atoms type A', '# atoms type B','# total atoms'})
legend boxoff


M = 5000; % total number of time steps
tspan = 0:(M-1);
y0=[N0 0 N0];
Time5=tic;
[t,y] = ode45(@(t,y)odefcn(y,lambda),tspan,y0);
Time6=toc(Time5);

f3 = figure;
figure(f3);
plot(t,y), box off
xlabel('Time (s)');
ylabel('Number of atoms');
title('Radioactive Decay (differential equation)','FontSize',14);
legend({'# atoms type A', '# atoms type B','# total atoms'})
legend boxoff

f4=figure;
figure(f4);
plot(t,sqrt((NC-NA).^2),'b',t,sqrt((y(:,1)-NA).^2),'red')
title('Variance from exact solution','FontSize',14);
xlabel('Time (s)');
ylabel('Number of atoms');
legend({'random function', 'differential equation'},'FontSize',14)


function dydt = odefcn(y, L)
dydt = zeros(3,1);
dydt(1) = -L*y(1);
dydt(2) = L*y(1);
dydt(3) = 0;
end

