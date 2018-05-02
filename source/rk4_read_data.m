
clc
clearvars
close all

%% Scalar system - C++/MATLAB

fileID1 = fileread('rk4_sys.txt') ;

data1 = textscan(fileID1, '%.5f\t%.5f', 'HeaderLines', 1);

t_vec1 = data1{:,1};
y_vec_C  = data1{:,2};

y0 = 1;

[t_vec1, y_vec23]  = ode23 (@fun_ode,t_vec1,y0);
[t_vec1, y_vec45]  = ode45 (@fun_ode,t_vec1,y0);
[t_vec1, y_vec113] = ode113(@fun_ode,t_vec1,y0);

figure('units','normalized','position',[.25 .25 .55 .4])
subplot(1,2,1);
plot(t_vec1,y_vec_C,'r','LineWidth',3);
hold on
plot(t_vec1,y_vec23,'b-.','LineWidth',1.5);
plot(t_vec1,y_vec45,'c-.','LineWidth',1.5);
plot(t_vec1,y_vec113,'m-.','LineWidth',1.5);
grid on
xlabel('t');
ylabel('y(t)');
legend({'rungekutta4','ode23','ode45','ode113'},'FontSize',11,'Location','best');
title('ODE solution RK4 and odexx');
%
subplot(1,2,2);
plot(t_vec1,y_vec_C - y_vec23,'b','LineWidth',1.5);
hold on
plot(t_vec1,y_vec_C - y_vec45,'c','LineWidth',1.5);
plot(t_vec1,y_vec_C - y_vec113,'m','LineWidth',1.5);
grid on
ylim([-10e-3 10e-3]);
xlabel('t');
ylabel('\Deltay = y_{C}(t) - y_{ML}(t)');
title('ODE \Delta, RK4 and odexx');
legend({'ML = ode23','ML = ode45','ML = ode113'},'FontSize',11,'Location','best');
saveas(gcf,'fig_rk4_sys.svg');

% some statistical data
vec_mean23  = mean(y_vec_C - y_vec23);
vec_mean45  = mean(y_vec_C - y_vec45);
vec_mean113 = mean(y_vec_C - y_vec113);
%
vec_std23  = std(y_vec_C - y_vec23);
vec_std45  = std(y_vec_C - y_vec45);
vec_std113 = std(y_vec_C - y_vec113);

%% Vector system - C++/MATLAB

fileID2 = fileread('rk4_sys_vec.txt') ;

DIM = 4;

data2 = textscan(fileID2, '%.5f\t%.5f\t%.5f\t%.5f\t%.5f', 'HeaderLines', 1);

t_vec2 = data2{:,1};
y_mat_C  = zeros(numel(t_vec2),DIM);
%
for j=1:DIM
    y_mat_C(:,j) = data2{:,j+1};
end

y0_vec = [0.3; 1.6; 0.9; 1.3];

[t_vec2, y_mat23]  = ode23(@(t,y) fun_odesys(t,y,DIM),t_vec2,y0_vec);
[t_vec2, y_mat45]  = ode45(@(t,y) fun_odesys(t,y,DIM),t_vec2,y0_vec);
[t_vec2, y_mat113] = ode113(@(t,y) fun_odesys(t,y,DIM),t_vec2,y0_vec);


for j=1:DIM
    
    figure('units','normalized','position',[.25 .25 .55 .4])
    subplot(1,2,1);
    plot(t_vec2,y_mat_C(:,j),'r','LineWidth',3);
    hold on
    plot(t_vec2,y_mat23(:,j),'b-.','LineWidth',1.5);
    plot(t_vec2,y_mat45(:,j),'c-.','LineWidth',1.5);
    plot(t_vec2,y_mat113(:,j),'m-.','LineWidth',1.5);
    grid on
    xlabel('t');
    ylabel(['y' num2str(j) ' (t)']);
    legend({'rungekutta4','ode23','ode45','ode113'},'FontSize',11,'Location','best');
    title(['ODE y' num2str(j) ' RK4 and odexx']);
    %
    subplot(1,2,2);
    plot(t_vec2,y_mat_C(:,j) - y_mat23(:,j),'b','LineWidth',1.5);
    hold on
    plot(t_vec2,y_mat_C(:,j) - y_mat45(:,j),'c','LineWidth',1.5);
    plot(t_vec2,y_mat_C(:,j) - y_mat113(:,j),'m','LineWidth',1.5);
    grid on
    xlabel('t');
    ylabel(['\Deltay = y' num2str(j) '_{C}(t) - y' num2str(j) '_{ML}(t)']);
    title(['ODE \Delta' num2str(j) ' , RK4 and odexx']);
    legend({'ML = ode23','ML = ode45','ML = ode113'},'FontSize',11,'Location','best');    
    %
    figname = 'fig_rk4_sys';
    fignumb = num2str(j);
    figexts = '.svg';
    saveas(gcf,strcat(figname,fignumb,figexts));
    
end

% Statistical data
mat_mean23  = mean(mean(y_mat_C - y_mat23));
mat_mean45  = mean(mean(y_mat_C - y_mat45));
mat_mean113 = mean(mean(y_mat_C - y_mat113));
%
mat_std23   = std(std(y_mat_C - y_mat23));
mat_std45   = std(std(y_mat_C - y_mat45));
mat_std113  = std(std(y_mat_C - y_mat113));

%% Function Declaration

function f = fun_ode(t,y)

    f = t*sin(y*t);
    
end

function f = fun_odesys(~,y,dim)

    f = zeros(dim,1);
    
    for j=0:dim-1
        f(j+1) = y(dim-j) - y(j+1)^2;
    end
    
end




