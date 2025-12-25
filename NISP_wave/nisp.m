

%%%%%%%%%% Code to Output MCS results from Python
clc
clear 



%%%%%%%%%%%%%% MCS Solutions from Python

u1m = load('./../outputs/3rv_sigma3_10000/u_MCS_mean.mat', 'u_MCS_mean');
u1s = load('./../outputs/3rv_sigma3_10000/u_MCS_sd.mat', 'u_MCS_sd');


u_MCS = u1m.u_MCS_mean;
u_MCS_sd = u1s.u_MCS_sd;


%%%%%%%% NISP solutions - Parallel %%%%%

u1m = load('./results/sigma3_d3l3/u_nisp_mean.mat', 'u_nisp_mean');
u1s = load('./results/sigma3_d3l3/u_nisp_sd.mat', 'u_nisp_sd');


u_nisp = u1m.u_nisp_mean;
u_nisp_sd = u1s.u_nisp_sd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% SSFEM Solver - Parallel %%%%%


filename_mean = fullfile('./../outputs/ssfem/sigma3_3rv_order4/mean_402.csv');
hpc_p = readtable(filename_mean);

filename_sd = fullfile('./../outputs/ssfem/sigma3_3rv_order4/sd_402.csv');
hpc_q = readtable(filename_sd);

t = hpc_p.Time;

t = t(:,1) + 1;

u_hpc = hpc_p.avg_u_;

u_hpc_sd = hpc_q.avg_u_;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




dt = 6.5e-3;

time = 0.01:dt:3+dt;

t_hpc = t.*dt;

nt = size(time,2);

pos = 403;


figure(24)
plot(time,u_MCS(pos,1:nt), time,u_nisp(pos,:),t_hpc, u_hpc, 'linewidth',2 )
xlabel('Time','fontSize',14);
ylabel('Pressure','fontSize',14);              
titlestring = ['MEAN - MCS and NISP'];
title(titlestring ,'fontsize',14);                            
h=gca; 
get(h,'FontSize');
set(h,'FontSize',14);
% fh = figure(length(plot_times)+1);
% set(fh, 'color', 'white'); 
xlim([0 3])
ylim([-inf inf])
legend('MCS-10000', 'NISP-L2-37','Intrusive-3rd order')


figure(25)
plot(time,u_MCS_sd(pos,1:nt),time,u_nisp_sd(pos,:),t_hpc, u_hpc_sd, 'linewidth',2)
xlabel('Time','fontSize',14);
ylabel('\sigma','fontSize',14);              
titlestring = ['SD - MCS and NISP'];
title(titlestring ,'fontsize',14);                            
h=gca; 
get(h,'FontSize');
set(h,'FontSize',14);
% fh = figure(length(plot_times)+1);
% set(fh, 'color', 'white'); 
xlim([0 3])
ylim([-inf inf])
legend('MCS-10000', 'NISP-L2-37','Intrusive-3rd order')



