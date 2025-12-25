
%%%%%%%% Code to calculate the confidence bound for MCS results %%%%%%%

dt = 6.5e-3;

samples = 10000;

y_mean = zeros(1,100);

n_exp = 10;

%%%%% K = 1,2,3 = 400 point - 150,300,462
%%%%% K = 4,5,6 = 401 point - 150,300,462
%%%%% K = 7,8,9 = 402 point - 150,300,462

k = 6;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
outputs_root = fullfile(repo_root,'outputs','3rv_sigma1');
filename = '3rv_pdf_%d'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding out the optimal x values 

path = fullfile(outputs_root,sprintf(filename,1));

file = load(fullfile(path,'u_pdf.mat'));

u_pdf = file.u_pdf;
   
[y,x] = ksdensity(u_pdf(k,:));

xmean = mean(x);
x_var = var(x);
 
sd_x = sqrt(x_var);
 
x_lim1 = xmean + 2*sd_x;
x_lim2 = xmean - 2*sd_x;

x = linspace(x_lim2,x_lim1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_exp

path = fullfile(outputs_root,sprintf(filename,i));

file = load(fullfile(path,'u_pdf.mat'));

u_pdf = file.u_pdf;


[y,x] = ksdensity(u_pdf(k,:), x);
 

y_mean = y_mean + y;

end

y_mean = y_mean/n_exp;


y_sd = zeros(1,100);


for i = 1:n_exp

path = fullfile(outputs_root,sprintf(filename,i));

file = load(fullfile(path,'u_pdf.mat'));

u_pdf = file.u_pdf;

[y,x] = ksdensity(u_pdf(k,:),x);


y_sd =  y_sd + (y - y_mean).^2;

end


y_sd = sqrt(y_sd/(n_exp -1));

% figure1 = figure;
% 
axes1 = axes('Parent',figure);
hold(axes1,'on');

YData1 = [y_mean+3*y_sd fliplr(y_mean-3*y_sd)];

XData1 = [x fliplr(x)];

% % Create patch
patch('Parent',axes1,'YData',YData1,'XData',XData1,...
    'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
    'EdgeColor','none');


% plot(x_mean,y_mean+y_sd,x_mean,y_mean-y_sd,'linewidth',2)
% patch([x_mean fliplr(x_mean)], [y_mean+2*y_sd fliplr(y_mean-2*y_sd)], 'y');


hold on

plot(x,y_mean,'k--', 'linewidth',2)
xlabel('x','fontSize',14);
ylabel('f(x)','fontSize',14);  
axis([-inf inf -inf inf])


    if(k == 1 || k == 4 || k == 7) 
        tt = 150*dt;
        titlestring = sprintf('PDF of solution process at time step %d ( T =%.3f s )',150, tt);
    elseif(k == 2 || k == 5 || k == 8)
        tt = 300*dt;
        titlestring = sprintf('PDF of solution process at time step %d ( T =%.3f s )',300, tt);
    else
        tt = 462*dt;
        titlestring = sprintf('PDF of solution process at time step %d ( T =%.3f s )',462, tt);
    end

title(titlestring ,'fontsize',14);                            
h=gca; 
get(h,'FontSize');
set(h,'FontSize',14);










