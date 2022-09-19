clear
clc 
set(0,'defaultfigurecolor','w');
%--------------------------------------------------------------------------
h0=load('参考线-h2dx.txt');%同性
h1=load('Irregular_var_dx_h2.txt');%同性

nstep = 6;
subplot(1,2,1)
%--------------------------------------------------------
  L=4; 
loglog(h1(1:nstep,1),abs(h1(1:nstep,L)),'ks-','LineWidth',1.2,'markersize',6)
hold on
loglog(h1(nstep+1:2*nstep,1),abs(h1(nstep+1:2*nstep,L)),'k^-','LineWidth',1.2,'markersize',5)
hold on
loglog(h1(2*nstep+1:3*nstep,1),abs(h1(2*nstep+1:3*nstep,L)),'ko-','LineWidth',1.2,'markersize',5)
hold on
%参考线///////////////
plot(h0(1:nstep,1),abs(h0(1:nstep,4)),'b-.','LineWidth',1.2,'markersize',6)
hold on
plot(h0(1:nstep,1),abs(h0(1:nstep,5)),'r--','LineWidth',1.2,'markersize',6)
hold on

ylabel('Error','FontSize',10,'FontName','Times New Roman');
xlabel('L/\Deltax','FontSize',10,'FontName','Times New Roman');
% legend('SPH \alpha_T/\alpha_L=0.1','FPM \alpha_T/\alpha_L=0.1');
legend('SPH','C-SPH','IC-SPH',' \propto\Deltax', ' \propto\Deltax^2','FontAngle','italic');
text(12,7e-4,'(a)', 'FontSize',12,'FontName','Times New Roman ');
grid on;
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','norm');
set(gcf,'Units','centimeters','Position',[1,1,10 10]);

% CPU--------------------------------------------------------------------------------------
subplot(1,2,2)
    L=4;
loglog(h1(1:nstep,7),abs(h1(1:nstep,L)),'ks-','LineWidth',1.2,'markersize',6)
hold on
loglog(h1(nstep+1:2*nstep,7),abs(h1(nstep+1:2*nstep,L)),'k^-','LineWidth',1.2,'markersize',5)
hold on
loglog(h1(2*nstep+1:3*nstep,7),abs(h1(2*nstep+1:3*nstep,L)),'ko-','LineWidth',1.2,'markersize',5)
hold on
%参考线///////////////
plot(h0(1:nstep,6),abs(h0(1:nstep,4)),'b-.','LineWidth',1.2,'markersize',6)
hold on
plot(h0(1:nstep,6),abs(h0(1:nstep,5)),'r--','LineWidth',1.2,'markersize',6)
hold on
ylabel('Error','FontSize',10,'FontName','Times New Roman');
xlabel('CPU Time (s)','FontSize',10,'FontName','Times New Roman');
legend('SPH','C-SPH','IC-SPH',' \propto\Deltax', ' \propto\Deltax^2','FontAngle','italic');
text(1.5,7e-4,'(b)', 'FontSize',12,'FontName','Times New Roman ');
set(gca, 'FontSize',12,'FontName','Times New Roman ','FontWeight','norm');
set(gcf,'Units','centimeters','Position',[1,1,22 10]); %19,14
grid on

