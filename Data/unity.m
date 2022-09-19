clc 
clear all
set(0,'defaultfigurecolor','w');
button=2; 
count=10; 
dt=3600;
maxtimestep=1000;
ntotal=160801; %160801,40401,10201,2601
p=6;
Dxx=1e-6;
Dyy=1e-6;
Dxy=0;
w=2;
x0=30;
y0=30;
c0=1;
C2=Dxx*Dyy-Dxy^2;
C3=Dxx+Dyy;

if (button==0) %//////////////////////////////////////////////////////////////////////////
elseif(button==1) %//////////////////////////////////////////////////////////////////////////////////
else %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ntotal=160401;  %251001，160801,40401,10201,2601
dt=3600;
results1=load('results_dx1.5_dx0.125_var0.mat');
results2=load('results_dx1.5_dx0.125_var0.25.mat');
results3=load('results_dx1.5_dx0.125_var0.5.mat');
results4=load('results_dx1.5_dx0.125_var1.mat');
results5=load('results_dx1.5_dx0.125_r.mat');

step=100; 
t=count*dt*step;
A1=2*t*Dyy+w^2; 
A2=2*t*Dxx+w^2;
A3=t*Dxy;
C2=Dxx*Dyy-Dxy^2;
C3=Dxx+Dyy;
C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);
%--------------------------------------------------------------------------
e=0.02;
yy=30;
k=0;
k2=0;
k3=0;
k4=0;
k5=0;
for n=1:ntotal
  
    if(  abs(results1.s(step,n,3)-yy) <= e  )  %选取1/2y处
        k=k+1;
        a(k,1)=results1.s(step,n,2); %x坐标数据
        a(k,2)=results1.s(step,n,3); %y坐标数据
        a(k,3)=c0*w^2/C4*exp((-( a(k,1)-x0)^2*A1-(a(k,2)-y0)^2*A2+4*( a(k,1)-x0)*( a(k,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
        
        a(k,4)=results1.s(step,n,6); 
        a(k,5)=results1.fc(step,n,6); 
        a(k,6)=results1.fcm(step,n,6); 
    end 
    %------------------------------------------------------
    if(  abs(results2.s(step,n,3)-yy) <= e  )  
         k2=k2+1;
         a2(k2,1)=results2.s(step,n,2); 
         a2(k2,2)=results2.s(step,n,3);
         a2(k2,3)=c0*w^2/C4*exp((-( a2(k2,1)-x0)^2*A1-(a2(k2,2)-y0)^2*A2+4*( a2(k2,1)-x0)*( a2(k2,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
        
         a2(k2,4)=results2.s(step,n,6);     
         a2(k2,5)=results2.fc(step,n,6); 
         a2(k2,6)=results2.fcm(step,n,6); 
    end
%     %-------------------------------------------------------
    if(  abs(results3.s(step,n,3)-yy) <= e  )  
         k3=k3+1;
         a3(k3,1)=results3.s(step,n,2); 
         a3(k3,2)=results3.s(step,n,3); 
         a3(k3,3)=c0*w^2/C4*exp((-( a3(k3,1)-x0)^2*A1-(a3(k3,2)-y0)^2*A2+4*( a3(k3,1)-x0)*( a3(k3,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
        
         a3(k3,4)=results3.s(step,n,6);  
         a3(k3,5)=results3.fc(step,n,6); 
         a3(k3,6)=results3.fcm(step,n,6); 
    end
    if(  abs(results4.s(step,n,3)-yy) <= e  ) 
         k4=k4+1;
         a4(k4,1)=results4.s(step,n,2); 
         a4(k4,2)=results4.s(step,n,3); 
         a4(k4,3)=c0*w^2/C4*exp((-( a4(k4,1)-x0)^2*A1-(a4(k4,2)-y0)^2*A2+4*( a4(k4,1)-x0)*( a4(k4,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
        
         a4(k4,4)=results4.s(step,n,6);  
         a4(k4,5)=results4.fc(step,n,6); 
         a4(k4,6)=results4.fcm(step,n,6); 
    end
       if(  abs(results5.s(step,n,3)-yy) <= e  )  %选取1/2y处
         k5=k5+1;
         a5(k5,1)=results5.s(step,n,2); 
         a5(k5,2)=results5.s(step,n,3); 
         a5(k5,3)=c0*w^2/C4*exp((-( a5(k5,1)-x0)^2*A1-(a5(k5,2)-y0)^2*A2+4*( a5(k5,1)-x0)*( a5(k5,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
        
         a5(k5,4)=results5.s(step,n,6);     
         a5(k5,5)=results5.fc(step,n,6);
         a5(k5,6)=results5.fdm(step,n,6); 
      end
end

% %//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
%四个拼图/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
b=sortrows(a,1); %根据x坐标排序
b2=sortrows(a2,1); 
b3=sortrows(a3,1);
b4=sortrows(a4,1);
b5=sortrows(a5,1);

c_max=0.6;
subplot(3,2,1)
xa=0:0.1:60;
ya=30;
cc=c0*w^2/C4.*exp((-( xa(:)-x0).^2*A1-(ya-y0)^2*A2+4*( xa(:)-x0).*( ya-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
p0 = plot(xa,cc,'k-','linewidth',1.5);
hold on
%------------------------------------------
%------------------------------------------
p1= plot(b(:,1),b(:,4),'bs','MarkerSize',3,'Markerfacecolor','b','linewidth',1.2) 
hold on
p2= plot(b(:,1),b(:,5),'^','color',[0,0.5,0],'MarkerSize',3,'Markerfacecolor','[0,0.5,0]')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
p3= plot(b(:,1),b(:,6),'ro','MarkerSize',3,'Markerfacecolor','r')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
xlim([15,45]);
ylim([0,c_max]);
ylabel('C/C_0','FontSize',10,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
lgd1=legend([p1,p2,p3,p0],'SPH','C-SPH','IC-SPH','Analytical');
set(lgd1,'box','off');
text(16,c_max-0.05,'(a) \sigma=0','FontSize',12,'FontName','Times New Roman','fontweight','bold');
set(gca, 'FontSize',10,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);

subplot(3,2,2)
xa=0:0.1:60;
ya=30;
cc=c0*w^2/C4.*exp((-( xa(:)-x0).^2*A1-(ya-y0)^2*A2+4*( xa(:)-x0).*( ya-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
p0=plot(xa,cc,'k-','linewidth',1.5);
hold on
%------------------------------------------
%------------------------------------------
maker_idx = 1:2:length(b2(:,1));
p1=plot(b2(:,1),b2(:,4),'bs','MarkerSize',3,'Markerfacecolor','b','linewidth',1.2,'MarkerIndices',maker_idx) 
hold on
p2=plot(b2(:,1),b2(:,5),'^','color',[0,0.5,0],'MarkerSize',3,'Markerfacecolor','[0,0.5,0]','MarkerIndices',maker_idx)  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
p3=plot(b2(:,1),b2(:,6),'ro','MarkerSize',3,'Markerfacecolor','r','MarkerIndices',maker_idx)  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
xlim([15,45]);
ylim([0,c_max]);
ylabel('C/C_0','FontSize',10,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
lgd1=legend([p1,p2,p3,p0],'SPH','C-SPH','IC-SPH','Analytical');
set(lgd1,'box','off');
text(16,c_max-0.05,'(b) \sigma=0.25\Deltax','FontSize',12,'FontName','Times New Roman','fontweight','bold');
set(gca, 'FontSize',10,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);

subplot(3,2,3)
% %解析解-------------------------------------
xa=0:0.1:60;
ya=30;
cc=c0*w^2/C4.*exp((-( xa(:)-x0).^2*A1-(ya-y0)^2*A2+4*( xa(:)-x0).*( ya-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
p0=plot(xa,cc,'k-','linewidth',1.5);
hold on
%------------------------------------------
p1=plot(b3(:,1),b3(:,4),'bs','MarkerSize',3,'Markerfacecolor','b','linewidth',1.2) 
hold on
p2=plot(b3(:,1),b3(:,5),'^','color',[0,0.5,0],'MarkerSize',3,'Markerfacecolor','[0,0.5,0]')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
p3=plot(b3(:,1),b3(:,6),'ro','MarkerSize',3,'Markerfacecolor','r')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
xlim([15,45]);
ylim([0,c_max]);
ylabel('C/C_0','FontSize',10,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
lgd1=legend([p1,p2,p3,p0],'SPH','C-SPH','IC-SPH','Analytical');
set(lgd1,'box','off');

text(16,c_max-0.05,'(c) \sigma=0.5\Deltax','FontSize',12,'FontName','Times New Roman','fontweight','bold');
set(gca, 'FontSize',10,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);

subplot(3,2,4)
%解析解-------------------------------------
xa=0:0.1:60;
ya=30;
cc=c0*w^2/C4.*exp((-( xa(:)-x0).^2*A1-(ya-y0)^2*A2+4*( xa(:)-x0).*( ya-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
p0=plot(xa,cc,'k-','linewidth',1.5);
hold on
%------------------------------------------
p1=plot(b4(:,1),b4(:,4),'bs','MarkerSize',3,'Markerfacecolor','b','linewidth',1.2) 
hold on
p2=plot(b4(:,1),b4(:,5),'^','color',[0,0.5,0],'MarkerSize',3,'Markerfacecolor','[0,0.5,0]')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
p3=plot(b4(:,1),b4(:,6),'ro','MarkerSize',3,'Markerfacecolor','r')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
xlim([15,45]);
ylim([0,c_max]);
ylabel('C/C_0','FontSize',10,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
lgd1=legend([p1,p2,p3,p0],'SPH','C-SPH','IC-SPH','Analytical');
set(lgd1,'box','off');
text(16,c_max-0.05,'(d) \sigma=\Deltax','FontSize',12,'FontName','Times New Roman','fontweight','bold');
set(gca, 'FontSize',10,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);

subplot(3,2,5.5)
% %解析解-------------------------------------
xa=0:0.1:60;
ya=30;
cc=c0*w^2/C4.*exp((-( xa(:)-x0).^2*A1-(ya-y0)^2*A2+4*( xa(:)-x0).*( ya-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%解析解
p0=plot(xa,cc,'k-','linewidth',1.5);
hold on
%------------------------------------------
p1=plot(b5(:,1),b5(:,4),'bs','MarkerSize',3,'Markerfacecolor','b','linewidth',1.2) 
hold on
p2=plot(b5(:,1),b5(:,5),'^','color',[0,0.5,0],'MarkerSize',3,'Markerfacecolor','[0,0.5,0]')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
p3=plot(b5(:,1),b5(:,6),'ro','MarkerSize',3,'Markerfacecolor','r')  %,'Color',[0,0.5,0],'color',[119,172,48]/255绿色
hold on
xlim([15,45]);
ylim([0,c_max]);
ylabel('C/C_0','FontSize',10,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
% lgd1=legend([p1,p2,p3,p0],'SPH','C-SPH','IC-SPH','Analytical');
% set(lgd1,'box','off');
text(16,c_max-0.05,'(e) Random','FontSize',12,'FontName','Times New Roman','fontweight','bold');
set(gca, 'FontSize',10,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);

set(gcf,'Units','centimeters','Position',[1 1 22 16]); %19,14
% set(gcf,'Units','centimeters','Position',[1 1 22 24]); %19,14
disp('done')
end  %///////////////////////////////////////////////////////////////////////////////////////////////////////////////////