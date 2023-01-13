clear, 
close all, 
N0 = 3e10;  % initial number of atoms

T1=7.6; %minute; Bi215
T2=36.1; %minute; Pb211
T3=2.14; %minute; Bi211
T4=4.77; %minute; Tl207
P3= 0.997; %probability Pb207

M = 200; % total number of time steps
tspan = 0:(M-1); 
y0=[N0 0 0 0 0];
[t,y] = ode45(@(t,y)ode_func(y,T1, T2, T3, T4, P3),tspan,y0);

for i = 1:M
    plot(t(1:i),y(1:i,:),'LineWidth',2)
    xlim([0 M])
    ylim([0 N0])
    xlabel('Time (min)');
    ylabel('Number of atoms');
    title('Radioactive Decay','FontSize',14);
    legend({'# atoms Bi215', '# atoms Pb211','# atoms Bi211', '# atoms Tl207','# atoms Pb207'},'FontSize',16),legend boxoff
    drawnow;
end



function dydt = ode_func(y, T1, T2, T3, T4,P3)
    dydt = zeros(5,1);
    dydt(1) = -y(1)/T1;
    dydt(2) = y(1)/T1 - y(2)/T2;
    dydt(3) = y(2)/T2 -y(3)/T3;
    dydt(4) = P3*y(3)/T3-y(4)/T4;
    dydt(5) = y(4)/T4+(1-P3)*y(3)/T3;
end
