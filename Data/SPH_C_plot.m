clear all
clc
set(0,'defaultfigurecolor','w');
t=3600*100*10; %3600*6*100;
Dxx=1e-6;
Dyy=1e-6;
Dxy=0;
w=2;
x0=30;
y0=30;
c0=1;
A1=2*t*Dyy+w^2;
A2=2*t*Dxx+w^2;
A3=t*Dxy;
C2=Dxx*Dyy-Dxy^2;
C3=Dxx+Dyy;
C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);

%2Dconcentration---------------------------------------------------------------
results1=load('results_dx1.5_dx0.125_var1.mat');
type=2; %1for C, 2for dC
step=100;
target(:,:)=results1.fdm(step,:,:);
cc=c0*w^2/C4*exp((-(target(:,2)-x0).^2*A1 -(target(:,3)-y0).^2*A2  +4.*(target(:,2)-x0).*(target(:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));

xv = linspace(5, 55, 1000);
yv = linspace(5, 55, 1000);
[Xm,Ym] = ndgrid(xv, yv);
x = target(:,2);
y = target(:,3);
z0=cc(:); %½âÎö½â

if(type==1) %c
z1 = results1.s(step,:,6);
z2=  results1.fc(step,:,6);
z3=  results1.fcm(step,:,6);
else %dc
z1 = results1.s(step,:,6)'-z0;
z2=  results1.fc(step,:,6)'-z0;
z3=  results1.fcm(step,:,6)'-z0;
end

max0=max(z0);
min0=min(z0);
min1 = min(z1);
min2 = min(z2);
min3 = min(z3);
max1 = max(z1);
max2 = max(z2);
max3 = max(z3);

nx=2;
ny=2;
fig_margin=0.08;
left_margin=0.1;
right_margin=0.1;
top_margin=0.05;
btm_margin=0.1;
fig_dx=(1-left_margin-right_margin-(nx-1)*fig_margin)/nx;
fig_dy=(1-top_margin-btm_margin-(ny-1)*fig_margin)/ny;

if(type==1)
subplot(2,2,1) %-------------------------
Zm = griddata(x, y, z0, Xm, Ym);
contourf(Xm, Ym, Zm, 50, 'LineColor','none');
% scatter(target(:,2),target(:,3),5,z0);
hold on
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
xlim([5,55]);
ylim([5,55]);
colormap jet
caxis([0 0.5]);
text(7,52,'(a1) Analytical','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin,btm_margin+(fig_dy+fig_margin),fig_dx,fig_dy]);

subplot(2,2,2) %---------------------------
Zm = griddata(x, y, z1, Xm, Ym);
contourf(Xm, Ym, Zm, 50,'LineColor','none');  % contourf(Xm, Ym, Zm, 5,'showtext','off', 'LineStyle', '--')
% scatter(target(:,2),target(:,3),5,z1);
hold on
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
xlim([5,55]);
ylim([5,55]);
colormap jet
caxis([0 0.5]);
% hcb=colorbar;
text(7,52,'(a2) SPH','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin+(fig_dx+fig_margin),btm_margin+(fig_dy+fig_margin),fig_dx,fig_dy]);

subplot(2,2,3) %---------------------------
Zm = griddata(x, y, z2, Xm, Ym);
contourf(Xm, Ym, Zm,50,'LineColor','none');
% scatter(target(:,2),target(:,3),5,z2);
hold on
xlim([5,55]);
ylim([5,55]);
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
colormap jet
caxis([0 0.5]);
% hcb=colorbar;
% title(hcb,'C/C_0')
text(7,52,'(a3) C-SPH','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin,btm_margin,fig_dx,fig_dy]);

subplot(2,2,4) %------------------------
Zm = griddata(x, y, z3, Xm, Ym);
contourf(Xm, Ym, Zm, 50,'LineColor','none')
% scatter(target(:,2),target(:,3),5,z3);
hold on
xlim([5,55]);
ylim([5,55]);
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
colormap jet
caxis([0 0.5]);
text(7,52,'(a4) IC-SPH','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin+(fig_dx+fig_margin),btm_margin,fig_dx,fig_dy]);
% colorbar('location','Eastoutside')
% colorbar('position',[0.5 1 1 100]);

top_margin = 0.1; % top margin
btm_margin = 0.1; % bottom margin
left_margin = 0.03;% left margin
right_margin = 0.2;% right margin
fig_margin = 0.08; % margin beween figures(sub) 
axes('position', [1-right_margin, btm_margin, 0.2, 1-(top_margin+btm_margin)]);
axis off;
colorbar;
caxis([0 0.5]);
hcb=colorbar;
title(hcb,'C/C_0');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gcf,'Units','centimeters','Position',[1 1 20 17]);

else %dc
subplot(2,2,1) %-------------------------
Zm = griddata(x, y, z1, Xm, Ym);
contourf(Xm, Ym, Zm, 50, 'LineColor','none');
% scatter(target(:,2),target(:,3),5,z0);
hold on
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
xlim([5,55]);
ylim([5,55]);
colormap jet
% caxis([0 0.5]);
caxis([-0.02 0.08]);
text(7,52,'(b1) SPH','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin,btm_margin+(fig_dy+fig_margin),fig_dx,fig_dy]);

subplot(2,2,2) %---------------------------
Zm = griddata(x, y, z2, Xm, Ym);
contourf(Xm, Ym, Zm, 50,'LineColor','none'); 
% scatter(target(:,2),target(:,3),5,z1);
hold on
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
xlim([5,55]);
ylim([5,55]);
colormap jet
caxis([0 0.5]);
caxis([-0.02 0.08]);
% hcb=colorbar;
text(7,52,'(b2) C-SPH','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin+(fig_dx+fig_margin),btm_margin+(fig_dy+fig_margin),fig_dx,fig_dy]);

subplot(2,2,3.5) %---------------------------
Zm = griddata(x, y, z3, Xm, Ym);
contourf(Xm, Ym, Zm,50,'LineColor','none');
% scatter(target(:,2),target(:,3),5,z2);
hold on
xlim([5,55]);
ylim([5,55]);
ylabel('y (m)','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
colormap jet
caxis([-0.02 0.08]);
% caxis([0 0.5]);
% hcb=colorbar;
% title(hcb,'C/C_0')
text(7,52,'(b3) IC-SPH','FontName','Times New Roman','FontSize',12,'Color','w','FontWeight','Bold');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gca,'position',[left_margin+0.5*(fig_dx+fig_margin),btm_margin,fig_dx,fig_dy]);
% left_margin+(fig_dx+fig_margin)

top_margin = 0.1; % top margin
btm_margin = 0.1; % bottom margin
left_margin = 0.03;% left margin
right_margin = 0.2;% right margin
fig_margin = 0.08; % margin beween figures(sub) 
axes('position', [1-right_margin+0.01, btm_margin, 0.2, 1-(top_margin+btm_margin)]);
axis off;
colorbar;
caxis([-0.02 0.08]);
hcb=colorbar;
title(hcb,'(C_N - C_A)/C_0')
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','normal');
set(gcf,'Units','centimeters','Position',[1 1 20 17]);
end