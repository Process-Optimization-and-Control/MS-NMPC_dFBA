% Plot FBA solution space from csv files
%reading files
% vg_val = readmatrix('vg_val.csv');
% vo_val = readmatrix('vo_val.csv');
vkmax = readmatrix('vk_max.csv');
vkmin = readmatrix('vk_min.csv');
vmax = readmatrix('v_max_1.csv');
vmin = readmatrix('v_min_1.csv');
vFBA = readmatrix('v_max_1_OF.csv');
vkFBA = readmatrix('vk_max_OF.csv');

vo_val=[15 14 13 12 11 10 9 8 7 6 5 3 2 1 0] ;
vg_val = [10 9.5 9 8.5 8 7.5 7 6.5 6 5.5 5.2 5 4.5 4 3.5];

% Acetate

f4 = figure(4);
% subplot(1,3,1)
surf(vg_val,vo_val,vmax)
xlabel('glucose uptake')
ylabel('oxygen uptake')
zlabel('Max Acetate flux')
% zlim([0 5])
title("FVA")
% title("\mu optimal")
% fontsize(f4, 20, "points");

f5 = figure(5);
% subplot(1,3,2)
surf(vg_val,vo_val,vkmax)
xlabel('glucose uptake')
ylabel('oxygen uptake')
zlabel('Max Acetate flux')
% zlim([0 5])
title("W = 0")
% title("\mu optimal")
% fontsize(f4, 20, "points");

% subplot(1,3,3)
f6 = figure(6);
surf(vg_val,vo_val,vkmax-vmax)
xlabel('glucose uptake')
ylabel('oxygen uptake')
zlabel('Max Acetate flux')
% zlim([0 5])
title("Residual")
% title("\mu optimal")
% fontsize(f4, 20, "points");
% 
% residual = norm(vkmax-vmax,'fro')

sum(sqrt(mean((vkmax - vmax).^2)))

% f4 = figure(4);
% % subplot(1,3,1)
% surf(vg_val,vo_val,vmin)
% xlabel('glucose uptake')
% ylabel('oxygen uptake')
% zlabel('Min Acetate flux')
% % zlim([0 5])
% title("FVA")
% % title("\mu optimal")
% % fontsize(f4, 20, "points");
% 
% f5 = figure(5);
% % subplot(1,3,2)
% surf(vg_val,vo_val,vkmin)
% xlabel('glucose uptake')
% ylabel('oxygen uptake')
% zlabel('Min Acetate flux')
% % zlim([0 5])
% title("W = 0")
% % title("\mu optimal")
% % fontsize(f4, 20, "points");
% 
% % subplot(1,3,3)
% f6 = figure(6);
% surf(vg_val,vo_val,vkmin-vmin)
% xlabel('glucose uptake')
% ylabel('oxygen uptake')
% zlabel('Min Acetate flux')
% % zlim([0 5])
% title("Residual")
% % title("\mu optimal")
% % fontsize(f4, 20, "points");

% % residual = norm(vkmin-vkmin,'fro')
% sum(sqrt(mean((vkmin - vmin).^2)))
% 
% f4 = figure(4);
% % subplot(2,1,1)
% surf(vg_val,vo_val,vkmax-vkmin)
% xlabel('glucose uptake')
% ylabel('oxygen uptake')
% zlabel('Max Acetate flux - Min acetate flux')
% zlim([0 5])
% title("W = +/- 1e-4")
% % title("\mu optimal")
% fontsize(f4, 20, "points");


% Growth rate
 f15 = figure(15);
% subplot(1,3,1)
surf(vg_val,vo_val,vFBA)
xlabel('glucose uptake')
ylabel('oxygen uptake')
zlabel('Growth rate')
% zlim([0 1])
title("LP")
fontsize(f15, 20, "points")
% 
% 
% 
% subplot(1,3,2)
% surf(vg_val,vo_val,(vkFBA-vkmax)./vkFBA)
% xlabel('glucose uptake')
% ylabel('oxygen uptake')
% zlabel('Growth rate residual (%)')
% zlim([0 0.2])
% title("W = +/- 2.8e-3")
% fontsize(f5, 20, "points")
% 
% 
% subplot(1,3,3)
f16 = figure(16);
surf(vg_val,vo_val,vkFBA)
xlabel('glucose uptake')
ylabel('oxygen uptake')
zlabel('Growth rate')
% zlim([0 1])
title("KKT")
% fontsize(f5, 20, "points")
sum(sqrt(mean((vkFBA - vFBA).^2)))