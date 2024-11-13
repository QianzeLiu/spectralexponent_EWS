%Coupled sea-ice-ocean model
clear
n2=1;n3=0.3;k=0.303;tt=200;ti=1;delta=0.43;h=0.5;R0=-0.1;B=0.45;L=1.25;F=1/28;
N1=100;
dt=0.1;t=100; N=t/dt;
xt=dt:dt:t;
I=zeros(N1,N);I(:,1)=1.1;
sigma=0.1;
R=0:-(0+0.5)/(N-1):-0.5;
for n=1:N1
    for i=2:N
        I(n,i)=I(n,i-1)+dt*(delta*tanh(I(n,i-1)/h)+(R0*heaviside(I(n,i-1))-B)*I(n,i-1)+L-F-1+R(i))+dt*sigma*randn/sqrt(dt);
    end
end
figure
plot(xt,mean(I))
%%
%Coupled sea-ice-ocean model
clear
n2=1;n3=0.3;k=0.303;tt=200;ti=1;delta=0.43;h=0.5;R0=-0.1;B=0.45;L=1.25;F=1/28;
r=-0.35:0.01:0.1;

for i=1:25
    R=r(i);
    syms I;
    y=delta*tanh(I/h)+(R0*heaviside(I)-B)*I+L-F-1+R;
    x1(i)=vpasolve(y==0,I,[-1.5,-0.47]);
end
for i=1:18
    R=r(i+7);
    syms I;
    y=delta*tanh(I/h)+(R0*heaviside(I)-B)*I+L-F-1+R;
    x2(i)=vpasolve(y==0,I,[-0.55,0.6]);
end
for i=1:39
    R=r(i+7);
    syms I;
    y=delta*tanh(I/h)+(R0*heaviside(I)-B)*I+L-F-1+R;
    x3(i)=vpasolve(y==0,I,[0.25,1.35]);
end
x1=eval(x1);
x2=eval(x2);
x3=eval(x3);

figure
plot(r(8:end),x3,'k',r(8:25),x2,'k:.',r(1:25),x1,'r','LineWidth',3),xlim([-0.35 0]),ylim([-1.49,1.49])
xlabel('The control parameter  (M)','FontWeight','bold'),ylabel('Sea ice  (I)','FontWeight','bold'),title('Bifurcation digram','FontName','Arial','fontweight','bold','FontSize',19)
legend('Ice covered','Unstable','Ice free','Location','northwest','fontweight','bold','FontSize',12),text(-0.35,1.6,'(a)','fontweight','bold',FontSize=19)
legend('boxoff')
set(gca,'fontweight','bold','linewidth',3,FontSize=19)

%%

step=199;W=200;
nfft=W;  Fs=W;       
c=0.01;d=.1;
a=W*c/2;b=W*d/2;
beta_per=zeros();beta_welch=zeros();beta_Yule=zeros();beta_burg=zeros();
for n=1:N1
    for m=1:N-step                      
        Y=I(n,m:m+step);


        window=boxcar(length(Y));
        [pxx,~]=periodogram(Y,window,nfft,Fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_per(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2+eps);

        % welch法
        
        win = hann(W/2); % 使用汉宁窗函数
        [pxx, ~] = pwelch(Y, win, [], [], Fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_welch(n,m)  = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);
 
        %Yule-Walker
        order = 10; 
        [pxx, ~] = pyulear(Y, order);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_Yule(n,m)  = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

        %Burg
        order = 10;
        [pxx, ~] = pburg(Y, order);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_burg(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

    end

end

%%

mp_per=mean(beta_per,1);
mp_welch=mean(beta_welch,1);
mp_Yule=mean(beta_Yule,1);
mp_burg=mean(beta_burg,1);

xl=100;
tau_x=1:(xl/dt-W);
tau1=corr(mp_per(1:(xl/dt-W))',tau_x','type','kendall');tau2=corr(mp_welch(1:(xl/dt-W))',tau_x','type','kendall');
tau3=corr(mp_Yule(1:(xl/dt-W))',tau_x','type','kendall');tau4=corr(mp_burg(1:(xl/dt-W))',tau_x','type','kendall');
figure
subplot(211)
ax1 = gca;
ax1.Position(2)=ax1.Position(2)-0.05;
set(gca,'Position',ax1.Position)
set(gcf,'position',[10 1 1050 900])
set(gca,'fontweight','bold','linewidth',3,FontSize=19)
set(gca,'XAxisLocation','top','YAxisLocation','right','XDir','reverse','xlim',[-10 0],'XTicklabels',[eval(vpa(-0.5:0.5/10:0 , 2))],'YTick',[]);
xlimits = get(ax1,'XLim');
xinc = (xlimits(2)-xlimits(1))/10;
set(ax1,'XTick',[xlimits(1):xinc:xlimits(2)])

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','left');
for n=1:100
    plot(xt,I(n,:))
    hold on
end
hold on
plot(xt,mean(I),'k','linewidth',3),xlabel(ax1,'Control parameter - M','FontWeight','bold'),xlabel(ax2,' '),set(ax2 ,'xticklabels', [])
xline(xl,'--k','linewidth',3)
ylabel('Sea Ice - I(t)','FontWeight','bold') ,xlim([0 100]),ylim([-2 1.5])
text(-0,2.2,'(b)','FontWeight','bold',FontSize=24)

axis normal
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(ax2,'XTick',[0:10:100])

subplot(212)
 err_per=std(beta_per);
 err_welch=std(beta_welch);
 err_Yule=std(beta_Yule);
 err_burg=std(beta_burg);


plot(xt(W+1:xl/dt),mp_per(1:(xl/dt-W)),'k',xt(W+1:xl/dt),mp_welch(1:(xl/dt-W)),'r',...
    xt(W+1:xl/dt),mp_Yule(1:(xl/dt-W)),'b',xt(W+1:xl/dt),mp_burg(1:(xl/dt-W)),'g','linewidth',3),xlim([0 100])
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[mp_per(1:(xl/dt-W))-abs(err_per(1:(xl/dt-W))),fliplr(mp_per(1:(xl/dt-W))+abs(err_per(1:(xl/dt-W)) ))],'r');
pa.FaceColor = "[0.5 .5 .5]";
pa.FaceAlpha = 0.2; 
pa.LineStyle = 'none'; 

pa1 = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[mp_welch(1:(xl/dt-W))-abs(err_welch(1:(xl/dt-W))),fliplr(mp_welch(1:(xl/dt-W))+abs(err_welch(1:(xl/dt-W)) ))],'r');
pa1.FaceColor = "#f06464";
pa1.FaceAlpha = 0.5;  
pa1.LineStyle = 'none';  

pa2 = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[mp_Yule(1:(xl/dt-W))-abs(err_Yule(1:(xl/dt-W))),fliplr(mp_Yule(1:(xl/dt-W))+abs(err_Yule(1:(xl/dt-W)) ))],'r');
pa2.FaceColor = "#2aaaef"; 
pa2.FaceAlpha = 0.5;  
pa2.LineStyle = 'none';  

pa3 = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[mp_burg(1:(xl/dt-W))-abs(err_burg(1:(xl/dt-W))),fliplr(mp_burg(1:(xl/dt-W))+abs(err_burg(1:(xl/dt-W)) ))],'r');
pa3.FaceColor = "#71c16f"; % 设置颜色
pa3.FaceAlpha = 0.3;  % 设置颜色透明度
 pa3.LineStyle = 'none';  % 设置误差带边界线为无



xlabel('Time','FontWeight','bold')
xline(xl,'--k','linewidth',3),
xline(xt(W+1),'--','color',[0.5 0.5 0.5],'linewidth',3)
text(1,min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)]),'- sliding window -',FontSize=14),text(-0,max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] ),'(c)','FontWeight','bold',FontSize=24)
ylim([min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)])-0.5 max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] )+0.7])
ylabel('Spectral exponent','FontWeight','bold') 
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
legend(['Periodogram      \tau = ',num2str(tau1)], ...
    ['Welch                 \tau = ',num2str(tau2)]', ...
    ['Yule-Walker       \tau = ',num2str(tau3)], ...
    ['Burg                   \tau = ',num2str(tau4)] , ...
    FontSize=12,Location='northeast');
legend('boxoff')
axis normal
set(gca,'XTick',[0:10:100])

%%
for n=1:N1
    tau_per(n)=corr(beta_per(n,1:(xl/dt-W))',tau_x','type','kendall');tau_wel(n)=corr(beta_welch(n,1:(xl/dt-W))',tau_x','type','kendall');
    tau_yul(n)=corr(beta_Yule(n,1:(xl/dt-W))',tau_x','type','kendall');tau_bur(n)=corr(beta_burg(n,1:(xl/dt-W))',tau_x','type','kendall');
end

tau=abs([tau_per;tau_wel;tau_yul;tau_bur]);
threshold=0.01:0.01:1;
color=['k','r','b','g'];
num=tau;
figure
set(gcf,'position',[10 100 800 800])
for n=1:4
    for i=1:length(threshold)
        th=threshold(i);
        num(n,i)=length(find(tau(n,:)>=th));
    end

    plot(threshold,num(n,:)/100,color(n),'linewidth',3)
    hold on
end
auc=sum(0.01.*num/100,2);
legend(['Periodogram  AUC = ',num2str(auc(1))], ...
    ['Welch             AUC = ',num2str(auc(2))], ...
    ['Yule-Walker   AUC = ',num2str(auc(3))],...
    ['Burg               AUC = ',num2str(auc(4))],FontSize=12,location='northeast')
legend('boxoff')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlabel('Kendall \tau threshold','FontWeight','bold')
ylabel('Proportion of regime shift detected','FontWeight','bold')
axis tight
text(-0,1.05,'(d)','FontWeight','bold',FontSize=24)


