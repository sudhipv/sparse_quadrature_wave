

%%%%%%% NISP PDF generation %%%%%%% 



%%%%% K = 1,2,3 = 400 point - 150,300,462
%%%%% K = 4,5,6 = 401 point - 150,300,462
%%%%% K = 7,8,9 = 402 point - 150,300,462


k = 5;
dt = 6.5e-3;
path = './results/3rv_d3l3_order4_pdf';
cd(path)

lmat = load('nisp_pdf_1.mat');

u_pdf = lmat.nisp_pdf_1;


[f_IN,xi_IN] = ksdensity(u_pdf(k,:));


figure(14)
plot(xi_IN,f_IN,'linewidth',2)
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

