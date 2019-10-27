clear 
clc
close all
tic
R=3390;
G=6.67259e-11;
g_mars_in=zeros(1,1);%Ԥ����neicun
g_mars_out=zeros(1,1);
V_mars_in=zeros(1,1);
V_mars_out=zeros(1,1);
M=zeros(1,1);
Rho=[7.266;7.266;7.036;7.036;4.220;4.220;4.007;4.007;3.599;3.599;2.850;2.850 ];
Rho=Rho.*1000;
r=[0;500;500;1468;1468;2033;2033;2360;2360;3280;3280;3390];
M(1,1)=(4/3).*pi.*Rho(1,1).*r(1,1)^3; %����������ʼֵ
for i=2:length(r) %�������ÿ���뾶ֵ����Ӧ������
    M(i,1)=M(i-1,1)+(4/3).*pi.*Rho(i,1).*(r(i,1)^3-r(i-1,1)^3);
end
%�ڲ���������ʼֵ
 g_mars_in(1,1)=(4/3).*pi.*G.*Rho(1,1).*r(1,1)*1000;
 
 %ѭ�����ڲ�������
for i=2:length(r) 
 g_mars_in(i,1)=G.*M(i,1).*1000/r(i,1).^2;
end
%����ڲ�������
V_out=zeros(1,1);
V_out(12,1)=0;
k=11;
while k>0  %�ⲿ���
    V_out(k,1)=10^(6).*G.*Rho(k,1).*(4.*pi/3).*(r(k+1,1)^2-r(k,1)^2)+V_out(k+1,1);
    k=k-1;
end
V_mars_in(1,1)=V_out(1,1);
for i=2:length(r) %���ڲ������ƿ��������֣��ⲿ���+�ڲ�����
     V_mars_in(i,1)=G.*M(i,1).*10^(6)/r(i,1)+V_out(i,1);
end
r_1=3390:3390*2;
%�ⲿ���������ⲿ������
for i=1:length(r_1)   
g_mars_out(i,1)=G.*M(12,1)*1000/r_1(1,i)^2;
V_mars_out(i,1)=G.*M(12,1)*10^(6)/r_1(1,i);
end

figure(1);
set(gcf,'unit','centimeters','position',[10.5,4.5,20,15.5]); %�Ի�ͼ����С����
[ax,h1,h2]=plotyy(r,Rho,r,g_mars_in);          %�Ȼ��ư뾶���ܶȹ�ϵͼ
set(gca,'box','off','Ytick',[]);        %�öδ��������Ϸ��̶��ص�������
set(ax(1),'ylim',[0,9000],'ytick',[0:1000:9000],'XAxisLocation','bottom'); %����ķ�Χ
set(ax(2),'ylim',[0,5],'ytick',[0:0.5:5],'XAxisLocation','bottom'); %����ķ�Χ 
view(90,-90);%��һ��������ʾ��z˳ʱ����תn�ȣ��ڶ���Ϊ��y
set(get(ax(1),'Ylabel'),'String','Density(kg/m^3)') %������������y��
set(get(ax(2),'Ylabel'),'String','Gravity of interior(m/s^2)') %�Ҳ�y��
set(gca,'Xgrid','on');%�԰뾶��Ϊ׼������
xlabel('Radius(km)');
set(gca,'XTick',[0:400:3500],'XLim',[0 3500]);
set(h1,'LineStyle','-','LineWidth',1.3)%��Ӧ��һ������
set(h2,'LineStyle','--','LineWidth',1.3)%��Ӧ�ڶ�������
lgd=legend('Rho(�����ܶ�)','g(�����ڲ�)','Location','west');%���ͼ������ͼ������������
lgd.FontSize = 12;                      %����ͼ�������С
legend('boxoff');                       %ȥ��ͼ����߿�
hold on
Ax1=axes('XAxisLocation','top','Color','none','YColor','none','XLim',[0 3500]); %��������ӻ����Ҳ�뾶������
xlabel(Ax1,'Radius(km)');
view(90,-90);%ԭʼֵΪview��0,90����ͼ����ת
set(gca,'XTick',[0:400:3500]); %�Ҳ�뾶��̶��޸�


figure(2);
set(gcf,'unit','centimeters','position',[10.5,4.5,20,15.5]);
[ax,h1,h2]=plotyy(r,V_mars_in,r,g_mars_in);          %�����ڲ�����λ��������ϵͼ
set(gca,'box','off','Ytick',[]); 
set(ax(1),'ylim',[0,3*10^(7)],'ytick',[0:5*10^(6):3*10^(7)]); %����ķ�Χ
set(ax(2),'ylim',[0,5],'ytick',[0:0.5:5]); %����ķ�Χ
view(90,-90);
set(get(ax(1),'Ylabel'),'String','Internal gravitational potential') %���y�� 
set(get(ax(2),'Ylabel'),'String','Gravity of interior(m/s^2)') %�Ҳ�y��
set(gca,'Xgrid','on');%�԰뾶��Ϊ׼������
set(gca,'XTick',[0:400:3500],'XLim',[0 3500]);
xlabel('Radius(km)');

set(h1,'LineStyle','-','LineWidth',1.3)%��Ӧ��һ�����ߵ�����y1
set(h2,'LineStyle','--','LineWidth',1.3)%��Ӧ�ڶ������ߵ�����y2 
lgd=legend('V(�����ڲ�)','g(�����ڲ�)','Location','west');
lgd.FontSize = 12;
legend('boxoff');
hold on
Ax1=axes('XAxisLocation','top','Color','none','YColor','none','XLim',[0 3500]); %��������ӻ����Ҳ�뾶������
xlabel(Ax1,'Radius(km)');
view(90,-90);   
set(gca,'XTick',[0:400:3500]);


figure(3);
set(gcf,'unit','centimeters','position',[10.5,4.5,20,15.5]);
[ax,h1,h2]=plotyy(r_1,V_mars_out,r_1,g_mars_out);          %�����ⲿ����λ��������ϵͼ
set(gca,'box','off','Ytick',[]); 
set(ax(1),'ylim',[0,3*10^(7)],'ytick',[0:5*10^(6):3*10^(7)]); %����ķ�Χ
set(ax(2),'ylim',[0,5],'ytick',[0:0.5:5]); %����ķ�Χ
view(90,-90);
set(get(ax(1),'Ylabel'),'String','External gravitational potential') %���y�� 
set(get(ax(2),'Ylabel'),'String','Gravity of External(m/s^2)') %�Ҳ�y��
set(gca,'Xgrid','on');%�԰뾶��Ϊ׼������
xlabel('Height(km)');
set(gca,'XTick',[3200:400:6800],'XLim',[3200 6800]);
set(h1,'LineStyle','-','LineWidth',1.3)%��Ӧ��һ�����ߵ�����y1
set(h2,'LineStyle','--','LineWidth',1.3)%��Ӧ�ڶ������ߵ�����y2 
lgd=legend('V(�����ⲿ)','g(�����ⲿ)','Location','best');
lgd.FontSize = 12;
legend('boxoff');
hold on
Ax1=axes('XAxisLocation','top','Color','none','YColor','none','XLim',[3200 3400*2]); %��������ӻ����Ҳ�뾶������
xlabel(Ax1,'Height(km)');
view(90,-90);
set(gca,'XTick',[3200:400:6800]);

toc


