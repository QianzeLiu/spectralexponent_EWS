clear
a=1.0;lamuda=0.12;p=1;Bc=10;B0=1;u=2;sigma1=0.01;sigma2=0.25;
dt=0.1;t=100; N=t/dt; N1=100;N2=40/dt;N3=70/dt;
R(1:N2)=1.7:(1.5-1.7)/(length(1:N2)-1):1.5;
R(N2:N3)=1.5:(1-1.5)/(length(N2:N3)-1):1;
R(N3:N)=1:(0.1-1)/(length(N3:N)-1):0.1;
w=zeros(N1,N);w(:,1)=1;
B=zeros(N1,N);B(:,1)=6.8;
xt=dt:dt:t;
for n=1:N1
    for i=2:N
        w(n,i)=w(n,i-1)+dt*(R(i-1)-a*w(n,i-1)-lamuda*w(n,i-1)*B(n,i-1))+dt*sigma1*randn/sqrt(dt);
        B(n,i)=B(n,i-1)+dt*(p*B(n,i-1)*(w(n,i-1)-B(n,i-1)/Bc)-u*B(n,i-1)/(B(n,i-1)+B0))-dt*sigma2*randn/sqrt(dt);
    end
end
figure
plot(xt,mean(B),'k')
xlabel('t'),ylabel(B(t))



%%
clear
clc
a=1.0;lamuda=0.12;p=1;Bc=10;B0=1;u=2;sigma1=0.01;sigma2=0.25;R=0.5:0.1:2.5;
A=zeros(3,length(R));C=A;
for i=1:length(R)
    syms w B

    f=R(i)-a*w-lamuda*w*B;
    g=p*B*(w-B/Bc)-u.*B/(B+B0);
    [x,y]=solve(f,g,w,B,'Real',true)
    %  [x,y]=solve(f,g,P,M,'MaxDegree',3)

    x=double(x(imag(double(x))==0));
    y=double(y(imag(double(y))==0));
    for n=1:length(x)
        A(n,i)=x(n);
    end
    for n=1:length(y)
        C(n,i)=y(n);
    end
end
double(x)
double(y)
A=sort(A,1,"descend");
C=sort(C,1,"descend");


figure
plot([1.05 R(7:end)],[C(1,7)/2+C(2,7)/2+0.1 C(1,7:end)],'k',[1.05 R(7:16)],[C(1,7)/2+C(2,7)/2-0.1 C(2,7:16) ],'k:.',R(1:16),C(3,1:16),'r','LineWidth',3)
ylim([-0.5,10.5])
xlabel('Rainfall rate (R)','FontWeight','bold'),ylabel('Vegetation biomass (B)','FontWeight','bold')
title('Bifurcation digram','FontName','Arial','FontSize',19,'fontweight','bold'),legend('Vegetated','Unstable','Bare','Location','northwest','FontSize',12)
text(.15,11.5,'(a)','fontweight','bold',FontSize=19)
set(gca,'fontweight','bold','linewidth',3,FontSize=19)
legend('boxoff')
figure
plot(R(1:16), A(1,1:16),'k',[1.05,R(7:16)], [A(2,7)/2+A(3,7)/2 A(2,7:15) 2 ],'k:.',[1.05 R(7:end)],[A(2,7)/2+A(3,7)/2-0.05 A(3,7:15) A(2,16) A(3,17:end)],'r','LineWidth',3)
ylim([0.4,2.1])
xlabel('Rainfall rate  (R)','FontWeight','bold'),ylabel(' Soil water  (w)','FontWeight','bold')
title('Bifurcation digram','FontName','Arial','FontSize',19,'fontweight','bold'),legend('Abundant','Unstable','Rare','Location','northwest','FontSize',12)
text(.15,2.2,'(b)','fontweight','bold',FontSize=19)
set(gca,'fontweight','bold','linewidth',3,FontSize=19)
legend('boxoff')


%%

step=199;W=200;
nfft=W; Fs=W;
c=0.01;d=0.1;
a=W*c/2;b=W*d/2;
beta_per=zeros();beta_welch=zeros();beta_Yule=zeros();beta_burg=zeros();
for n=1:N1
    for m=1:N-step
        Y=B(n,m:m+step);

        window=boxcar(length(Y));
        [pxx,~]=periodogram(Y,window,nfft,Fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_per(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2+eps);

        % welch法

        win = hann(W/2);
        [pxx, ~] = pwelch(Y, win, [], [], Fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_welch(n,m)  = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

        % Yule-Walker
        order = 10;
        [pxx, ~] = pyulear(Y, order);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_Yule(n,m)  = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

        % Burg
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

xl=80;

mp_per=mean(beta_per,1);

mp_welch=mean(beta_welch,1);

mp_Yule=mean(beta_Yule,1);

mp_burg=mean(beta_burg,1);


tau_x=1:(xl/dt-W);
tau1=corr(mp_per(1:(xl/dt-W))',tau_x','type','kendall');tau2=corr(mp_welch(1:(xl/dt-W))',tau_x','type','kendall');
tau3=corr(mp_Yule(1:(xl/dt-W))',tau_x','type','kendall');tau4=corr(mp_burg(1:(xl/dt-W))',tau_x','type','kendall');

figure
set(gcf,'position',[10 1 1050 900])
subplot(211)
ax1 = gca;
ax1.Position(2)=ax1.Position(2)-0.05;
set(gca,'Position',ax1.Position)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(gca,'XAxisLocation','top','YAxisLocation','right','XDir','reverse','xlim',[0 10],'XTicklabels',[0.1,0.4,0.7 1 1.16 1.33 1.5 1.55  1.6 1.65 1.7 ],'YTick',[]);
xlimits = get(ax1,'XLim');
xinc = (xlimits(2)-xlimits(1))/10;
set(ax1,'XTick',[xlimits(1):xinc:xlimits(2)])

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','left');
for n=1:100
    plot(xt,B(n,:))
    hold on
end
hold on
plot(xt,mean(B),'k','linewidth',3),xlabel(ax1,'Rainfall rate - R','FontWeight','bold'),xlabel(ax2,' '),set(ax2 ,'xticklabels', [])
xline(xl,'--k','linewidth',3)
ylabel('Biomass density - B(t)','FontWeight','bold') ,xlim([0 100]),ylim([-0.5 7.9]),
text(0,9.5,'(c)','FontWeight','bold',FontSize=24)
set(ax2,'XTick',[0:10:100])
axis normal
set(gca,'fontweight','bold','linewidth',4,FontSize=19)


subplot(212)
plot(xt(W+1:xl/dt),mp_per(1:(xl/dt-W)),'k',xt(W+1:xl/dt),mp_welch(1:(xl/dt-W)),'r',...
    xt(W+1:xl/dt),mp_Yule(1:(xl/dt-W)),'b',xt(W+1:xl/dt),mp_burg(1:(xl/dt-W)),'g','linewidth',3),xlim([0 100])

err_per=std(beta_per);
err_welch=std(beta_welch);
err_Yule=std(beta_Yule);
err_burg=std(beta_burg);

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
pa3.FaceColor = "#71c16f";
pa3.FaceAlpha = 0.3; 
pa3.LineStyle = 'none'; 



xlabel('Time','FontWeight','bold')
xline(xl,'--k','linewidth',3),
xline(xt(W+1),'--','color',[0.5 0.5 0.5],'linewidth',3)
text(1,min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)]),'- sliding window -',FontSize=14),text(-0,max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] ),'(d)','FontWeight','bold',FontSize=24)
ylim([min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)])-0.5 max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] )+1])
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
set(gcf,'position',[10 100 700 700])
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
    ['Burg               AUC = ',num2str(auc(4))],FontSize=12)
legend('boxoff')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlabel('Kendall \tau threshold','FontWeight','bold')
ylabel('Proportion of regime shift detected','FontWeight','bold')
axis tight
text(0,1.05,'(e)','FontWeight','bold',FontSize=24)


