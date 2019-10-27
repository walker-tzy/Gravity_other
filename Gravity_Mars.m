clear 
clc
close all
tic
R=3390;
G=6.67259e-11;
g_mars_in=zeros(1,1);%预分配neicun
g_mars_out=zeros(1,1);
V_mars_in=zeros(1,1);
V_mars_out=zeros(1,1);
M=zeros(1,1);
Rho=[7.266;7.266;7.036;7.036;4.220;4.220;4.007;4.007;3.599;3.599;2.850;2.850 ];
Rho=Rho.*1000;
r=[0;500;500;1468;1468;2033;2033;2360;2360;3280;3280;3390];
M(1,1)=(4/3).*pi.*Rho(1,1).*r(1,1)^3; %计算质量初始值
for i=2:length(r) %迭代求出每个半径值所对应的质量
    M(i,1)=M(i-1,1)+(4/3).*pi.*Rho(i,1).*(r(i,1)^3-r(i-1,1)^3);
end
%内部引力场初始值
 g_mars_in(1,1)=(4/3).*pi.*G.*Rho(1,1).*r(1,1)*1000;
 
 %循环求内部引力场
for i=2:length(r) 
 g_mars_in(i,1)=G.*M(i,1).*1000/r(i,1).^2;
end
%求解内部引力势
V_out=zeros(1,1);
V_out(12,1)=0;
k=11;
while k>0  %外部球壳
    V_out(k,1)=10^(6).*G.*Rho(k,1).*(4.*pi/3).*(r(k+1,1)^2-r(k,1)^2)+V_out(k+1,1);
    k=k-1;
end
V_mars_in(1,1)=V_out(1,1);
for i=2:length(r) %将内部引力势看成两部分，外部球壳+内部球体
     V_mars_in(i,1)=G.*M(i,1).*10^(6)/r(i,1)+V_out(i,1);
end
r_1=3390:3390*2;
%外部引力场和外部引力势
for i=1:length(r_1)   
g_mars_out(i,1)=G.*M(12,1)*1000/r_1(1,i)^2;
V_mars_out(i,1)=G.*M(12,1)*10^(6)/r_1(1,i);
end

figure(1);
set(gcf,'unit','centimeters','position',[10.5,4.5,20,15.5]); %对绘图区大小控制
[ax,h1,h2]=plotyy(r,Rho,r,g_mars_in);          %先绘制半径和密度关系图
set(gca,'box','off','Ytick',[]);        %该段代码解决了上方刻度重叠的问题
set(ax(1),'ylim',[0,9000],'ytick',[0:1000:9000],'XAxisLocation','bottom'); %右轴的范围
set(ax(2),'ylim',[0,5],'ytick',[0:0.5:5],'XAxisLocation','bottom'); %左轴的范围 
view(90,-90);%第一个参数表示沿z顺时针旋转n度，第二个为沿y
set(get(ax(1),'Ylabel'),'String','Density(kg/m^3)') %函数句柄，左侧y轴
set(get(ax(2),'Ylabel'),'String','Gravity of interior(m/s^2)') %右侧y轴
set(gca,'Xgrid','on');%以半径轴为准画网格
xlabel('Radius(km)');
set(gca,'XTick',[0:400:3500],'XLim',[0 3500]);
set(h1,'LineStyle','-','LineWidth',1.3)%对应第一条曲线
set(h2,'LineStyle','--','LineWidth',1.3)%对应第二条曲线
lgd=legend('Rho(火星密度)','g(火星内部)','Location','west');%添加图例，将图例放置在西侧
lgd.FontSize = 12;                      %设置图例字体大小
legend('boxoff');                       %去掉图例外边框
hold on
Ax1=axes('XAxisLocation','top','Color','none','YColor','none','XLim',[0 3500]); %用坐标叠加绘制右侧半径坐标轴
xlabel(Ax1,'Radius(km)');
view(90,-90);%原始值为view（0,90），图像旋转
set(gca,'XTick',[0:400:3500]); %右侧半径轴刻度修改


figure(2);
set(gcf,'unit','centimeters','position',[10.5,4.5,20,15.5]);
[ax,h1,h2]=plotyy(r,V_mars_in,r,g_mars_in);          %绘制内部引力位与引力关系图
set(gca,'box','off','Ytick',[]); 
set(ax(1),'ylim',[0,3*10^(7)],'ytick',[0:5*10^(6):3*10^(7)]); %右轴的范围
set(ax(2),'ylim',[0,5],'ytick',[0:0.5:5]); %左轴的范围
view(90,-90);
set(get(ax(1),'Ylabel'),'String','Internal gravitational potential') %左侧y轴 
set(get(ax(2),'Ylabel'),'String','Gravity of interior(m/s^2)') %右侧y轴
set(gca,'Xgrid','on');%以半径轴为准画网格
set(gca,'XTick',[0:400:3500],'XLim',[0 3500]);
xlabel('Radius(km)');

set(h1,'LineStyle','-','LineWidth',1.3)%对应第一条曲线的线性y1
set(h2,'LineStyle','--','LineWidth',1.3)%对应第二条曲线的线性y2 
lgd=legend('V(火星内部)','g(火星内部)','Location','west');
lgd.FontSize = 12;
legend('boxoff');
hold on
Ax1=axes('XAxisLocation','top','Color','none','YColor','none','XLim',[0 3500]); %用坐标叠加绘制右侧半径坐标轴
xlabel(Ax1,'Radius(km)');
view(90,-90);   
set(gca,'XTick',[0:400:3500]);


figure(3);
set(gcf,'unit','centimeters','position',[10.5,4.5,20,15.5]);
[ax,h1,h2]=plotyy(r_1,V_mars_out,r_1,g_mars_out);          %绘制外部引力位与引力关系图
set(gca,'box','off','Ytick',[]); 
set(ax(1),'ylim',[0,3*10^(7)],'ytick',[0:5*10^(6):3*10^(7)]); %右轴的范围
set(ax(2),'ylim',[0,5],'ytick',[0:0.5:5]); %左轴的范围
view(90,-90);
set(get(ax(1),'Ylabel'),'String','External gravitational potential') %左侧y轴 
set(get(ax(2),'Ylabel'),'String','Gravity of External(m/s^2)') %右侧y轴
set(gca,'Xgrid','on');%以半径轴为准画网格
xlabel('Height(km)');
set(gca,'XTick',[3200:400:6800],'XLim',[3200 6800]);
set(h1,'LineStyle','-','LineWidth',1.3)%对应第一条曲线的线性y1
set(h2,'LineStyle','--','LineWidth',1.3)%对应第二条曲线的线性y2 
lgd=legend('V(火星外部)','g(火星外部)','Location','best');
lgd.FontSize = 12;
legend('boxoff');
hold on
Ax1=axes('XAxisLocation','top','Color','none','YColor','none','XLim',[3200 3400*2]); %用坐标叠加绘制右侧半径坐标轴
xlabel(Ax1,'Height(km)');
view(90,-90);
set(gca,'XTick',[3200:400:6800]);

toc


