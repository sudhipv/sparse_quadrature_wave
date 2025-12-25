




%%%%%%%%%% Code to Output MCS results from Python
clc
clear 
% ExtractGmsh
file = fullfile('./../meshData/square_acoustic.msh');


%%%%% SSFEM WAVE SOLVER RESULTS - CSV FILE


%%%%%%%%%%%%%%%%%

ExtractGmsh(file);

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';

%%%%%%%%%%%%%% MCS Solutions from Python

% u1m = load('./membrane/3000/u_MCS_mean.mat', 'u_MCS_mean');
% u1s = load('./membrane/3000/u_MCS_sd.mat', 'u_MCS_sd');
% 
% 
% u_MCS_mean = u1m.u_MCS_mean;
% u_MCS_sd = sqrt(u1s.u_MCS_sd);



% u1m = load('./gaussian/3000/u_MCS_mean.mat', 'u_MCS_mean');
% u1s = load('./gaussian/3000/u_MCS_sd.mat', 'u_MCS_sd');
% 
% 
% u_MCS_mean_g = u1m.u_MCS_mean;
% u_MCS_sd_g = sqrt(u1s.u_MCS_sd);


% nS1 = load('./gaussian/nS.mat', 'nS');

% 
u1m10 = load('./gaussian/3s_3rv_check/u_MCS_mean.mat', 'u_MCS_mean');
u1s10 = load('./gaussian/3s_3rv_check/u_MCS_sd.mat', 'u_MCS_sd');


u_MCS_mean = u1m10.u_MCS_mean;
u_MCS_sd = u1s10.u_MCS_sd;


% u_MCS_sd = sqrt(u_MCS_sd);


% % % %%%%%%% 20000 sample
% u1m20 = load('./gaussian/3s_5rv_20000/u_MCS_mean.mat', 'u_MCS_mean');
% u1s20 = load('./gaussian/3s_5rv_20000/u_MCS_sd.mat', 'u_MCS_sd');
% 
% 
% u_MCS_mean_20 = u1m20.u_MCS_mean;
% u_MCS_sd_20 = u1s20.u_MCS_sd;


%%%%%%%%%%%%%%  SSFEM Solutions from HPC

% filename_mean = fullfile('./../ssfem_solver/13k_3rv_3s/mean_402.csv');
% hpc_p = readtable(filename_mean);
% 
% filename_sd = fullfile('./../ssfem_solver/13k_3rv_3s/sd_402.csv');
% hpc_q = readtable(filename_sd);

% t = hpc_p.Time;
% 
% t = t(:,1) + 1;
% 
% u_hpc = hpc_p.avg_u_;
% 
% u_hpc_sd = hpc_q.avg_u_;
% 



filename_mean = fullfile('./../ssfem_solver/13k_3rv_3s/mean_400.csv');
hpc_p = readtable(filename_mean);

filename_sd = fullfile('./../ssfem_solver/13k_3rv_3s/sd_400.csv');
hpc_q = readtable(filename_sd);

t = hpc_p.Time;

t = t(:,1) + 1;

u_hpc = hpc_p.avg_u_;

u_hpc_sd = hpc_q.avg_u_;




% filename_mean = fullfile('./../ssfem_solver/37k_3rv_3s/mean_668.csv');
% hpc_p = readtable(filename_mean);
% 
% filename_sd = fullfile('./../ssfem_solver/37k_3rv_3s/sd_668.csv');
% hpc_q = readtable(filename_sd);
% 
% t_37k = hpc_p.Time;
% 
% t_37k = t_37k(:,1) + 1;
% 
% u_hpc_37k = hpc_p.avg_u_;
% 
% u_hpc_sd_37k = hpc_q.avg_u_;







% filename_mean = fullfile('./../ssfem_solver/13k_damped_CFL65/mean_401.csv');
% hpc_p = readtable(filename_mean);
% 
% filename_sd = fullfile('./../ssfem_solver/13k_damped_CFL65/sd_401.csv');
% hpc_q = readtable(filename_sd);
% 
% 
% t = hpc_p.Time;
% 
% t = t(:,1) + 1;
% 
% u_hpc = hpc_p.avg_u_;
% 
% u_hpc_sd = hpc_q.avg_u_;
% 
% 
% 
% filename_mean = fullfile('./../ssfem_solver/13k_damped_CFL65_4thorder/mean_401.csv');
% hpc_p4th = readtable(filename_mean);
% 
% filename_sd = fullfile('./../ssfem_solver/13k_damped_CFL65_4thorder/sd_401.csv');
% hpc_4th = readtable(filename_sd);
% 
% 
% u_hpc_4th = hpc_p4th.avg_u_;
% u_hpc_sd4th = hpc_4th.avg_u_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% nt = size(u_MCS_mean_membrane,2);

dt = 6.5e-3;

time = dt:dt:3+dt;

t_hpc = t.*dt;

% t_hpc_37k = t_37k.*3.9e-3;

pos = 401;

% Mean Solution
% for j =1:nt
%     
%     
%      s1 = pdesurf(p,t,u_MCS_mean_membrane(:,j));
%      axis([0 1 0 1 -1 1]);
%      zlim([-1 1])
%      hold on
%      colorbar
%      pause(0.05);
%      
%          
%     if(j<nt)
%          delete(s1)
%     end
% 
% end
% 
% 
% for k =1:nt
%     
%     
%      s2 = pdesurf(p,t,u_MCS_sd_membrane(:,k));
%      axis([0 1 0 1 -1 1]);
%      zlim([-1 1])
%      hold on
%      colorbar
%      pause(0.05);
%      
%          
%     if(k<nt)
%          delete(s2)
%     end
% 
% end
% time, u_MCS_mean_20(pos,:)

% time, u_MCS_sd_20(pos,:)
% 
% t_hpc,u_hpc,    t_hpc,u_hpc_sd

% time, u_MCS_mean(pos,:),

figure(24)
plot(time, u_MCS_mean(pos, :),t_hpc,u_hpc, 'linewidth',2 )
xlabel('Time','fontSize',14);
ylabel('Pressure','fontSize',14);  
axis([-inf inf -inf inf])
titlestring = ['MEAN - MCS'];
title(titlestring ,'fontsize',14);                            
h=gca; 
get(h,'FontSize');
set(h,'FontSize',14);
% fh = figure(length(plot_times)+1);
% set(fh, 'color', 'white'); 
legend

% time, u_MCS_sd(pos,:)
% time, u_MCS_sd(pos,:)

figure(25)
plot(time, u_MCS_sd(pos, :),t_hpc,u_hpc_sd, 'linewidth',2)
xlabel('Time (s)','fontSize',14);
ylabel('Standard Deviation (\sigma)','fontSize',14);              
% titlestring = ['SD - MCS and SSFEM SOLVER'];
% title(titlestring ,'fontsize',14); 
axis([-inf inf -inf inf])
h=gca; 
get(h,'FontSize');
set(h,'FontSize',14);
% fh = figure(length(plot_times)+1);
% set(fh, 'color', 'white'); 

legend
% legend ('MCS-10000', 'MCS-20000', 'SSFEM - 3rd Order')



%%%%% CoV %%%%%%%%%%%%%

% u_hpc_cov = u_hpc_sd./u_hpc;
% 
% figure(26)
% plot(t_hpc,u_hpc_cov, 'linewidth',2)
% xlabel('Time','fontSize',14);
% ylabel('Pressure','fontSize',14);              
% titlestring = ['CoV for SSFEM Solver'];
% title(titlestring ,'fontsize',14);                            
% h=gca; 
% get(h,'FontSize');
% set(h,'FontSize',14);
% % fh = figure(length(plot_times)+1);
% % set(fh, 'color', 'white'); 
% legend ('SSFEM-HPC')
% 
% 




%%%%%%%%Error plot %%%%%%%


% % % error_mean = (u_MCS_mean(pos, :) - u_hpc');
% % % 
% % % error_mean = error_mean./u_MCS_mean(pos, :);
% % % 
% % % error_sd = (u_MCS_sd(pos,:) - u_hpc_sd');
% % % 
% % % error_sd = error_sd./u_MCS_sd(pos,:);
% % % 
% % % figure(26)
% % % plot(t_hpc,error_mean, 'linewidth',2)
% % % xlabel('Time','fontSize',14);
% % % ylabel('Error','fontSize',14);              
% % % % titlestring = ['CoV for SSFEM Solver'];
% % % % title(titlestring ,'fontsize',14);                            
% % % h=gca; 
% % % get(h,'FontSize');
% % % set(h,'FontSize',14);
% % % % fh = figure(length(plot_times)+1);
% % % % set(fh, 'color', 'white'); 
% % % legend ('Error in SSFEM-HPC')
% % % 
% % % 
% % % 
% % % figure(27)
% % % plot(t_hpc,error_sd, 'linewidth',2)
% % % xlabel('Time','fontSize',14);
% % % ylabel('Error','fontSize',14);              
% % % % titlestring = ['CoV for SSFEM Solver'];
% % % % title(titlestring ,'fontsize',14);                            
% % % h=gca; 
% % % get(h,'FontSize');
% % % set(h,'FontSize',14);
% % % % fh = figure(length(plot_times)+1);
% % % % set(fh, 'color', 'white'); 
% % % legend ('Error in SSFEM-HPC')







