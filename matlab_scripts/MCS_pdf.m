


%%%%%% Code to extract the PDF from sample points in MCS


script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
path = fullfile(repo_root,'outputs','3rv_sigma3','3rv_sigma3_pdf_20');
% 
% %%%%%%%%%%%%%%%%%%% Preprocess %%%%%%%%%%%%%%%%%%%%%%%%
% 
N = 400;
samples = 10000;

chunks = samples/N;


u_pdf = zeros(9,samples);

for i = 1:N

filename = fullfile(path,sprintf('u_pdf_%d.mat',i-1));

file = load(filename);

u_pdf(:,(i-1)*chunks+1 :  i*chunks) = file.u_pdf;

i = 1+1;

delete(filename)

end

save(fullfile(path,'u_pdf.mat'),'u_pdf')





%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%

% % % dt = 6.5e-3;
% % % 
% % % cd(path)
% % % 
% % % file = load('u_pdf.mat');
% % % 
% % % u_pdf = file.u_pdf;
% % % 
% % % 
% % % for k = 1:9
% % %  
% % %     figure(k)
% % %     
% % %     [y,x] = ksdensity(u_pdf(k,:));
% % %     
% % %     plot(x,y,'linewidth',2)
% % %     xlabel('x','fontSize',14);
% % %     ylabel('f(x)','fontSize',14);  
% % %     axis([-inf inf -inf inf])
% % %     if(k == 1 || k == 4 || k == 7) 
% % %         tt = 150*dt;
% % %         titlestring = sprintf('PDF of solution process at time step %d ( T =%.3f s )',150, tt);
% % %     elseif(k == 2 || k == 5 || k == 8)
% % %         tt = 300*dt;
% % %         titlestring = sprintf('PDF of solution process at time step %d ( T =%.3f s )',300, tt);
% % %     else
% % %         tt = 462*dt;
% % %         titlestring = sprintf('PDF of solution process at time step %d ( T =%.3f s )',462, tt);
% % %     end
% % %     title(titlestring ,'fontsize',14);                            
% % %     h=gca; 
% % %     get(h,'FontSize');
% % %     set(h,'FontSize',14);
% % %     
% % % end
% % % 



