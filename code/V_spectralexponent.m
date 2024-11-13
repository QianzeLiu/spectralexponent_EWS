
clear
r=1;K=10;V0=1;sigma=0.25;
dt=0.1;t=100; N=t/dt; N1=100;
V=zeros(N1,N);  V(:,1)=8.89;
xt=dt:dt:t;

C=1:(3-1)/(N-1):3;
for n=1:N1
    for i=2:N
        V(n,i)=V(n,i-1)+dt*(r*V(n,i-1)*(1-V(n,i-1)/K)-C(i-1)*V(n,i-1)*V(n,i-1)/(V(n,i-1)*V(n,i-1)+V0*V0))+dt*sigma*randn/sqrt(dt);
    end
end

figure
plot(xt,mean(V),'k'),xlabel('Time'),ylabel('V')




%%

% dV/dt=r*V*(1-V/K)-cV.*V/(V.*V+V0*V0)
clear
r=1;K=10;V0=1;c=1.5:0.05:3;
for i=1:23
    C=c(i);
    syms x;
    y=r*x.*(1-x./K)-C*x.^2/(x.^2+V0*V0);
    x1(i)=vpasolve(y==0,x,[5,8.2]);
end
for i=7:23
    C=c(i);
    syms x;
    y=r*x.*(1-x./K)-C*x.^2/(x.^2+V0*V0);
    x2(i-6)=vpasolve(y==0,x,[1.28,4.9]);
end
for i=7:31
    C=c(i);
    syms x;
    y=r*x.*(1-x./K)-C*x.^2/(x.^2+V0*V0);
    x3(i-6)=vpasolve(y==0,x,[0.35,1.1]);
end
x1=eval(x1);
x2=eval(x2);
x3=eval(x3);
figure
set(gcf,'position',[10 10 500 500])
plot(c(1:23),x1,'k',c(7:23),x2,'k:.',c(7:end),x3,'r','LineWidth',3)
ylim([0,9])
xlabel('Maximum grazing rate  (c)','FontWeight','bold'),ylabel('Vegetation biomass  (V)','FontWeight','bold'),
title('Bifurcation digram','FontName','Arial','FontSize',19,'fontweight','bold')
legend('Vegetated','Unstable','Bare','Location','northeast','FontSize',12)
legend('boxoff')
text(1.3,9.5,'(c)','fontweight','bold',FontSize=19)
set(gca,'fontweight','bold','linewidth',3,FontSize=19)
%%
step=199;W=200;
nfft=W;      
c=0.01;d=0.1;
a=W*c/2;b=W*d/2;
Fs=W;
beta_per=zeros();beta_welch=zeros();beta_Yule=zeros();beta_burg=zeros();

for n=1:N1
    for m=1:N-step                     
        Y=V(n,m:m+step);

        %periodogram
        window=boxcar(length(Y));
        [pxx,~]=periodogram(Y,window,nfft,Fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_per(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2+eps);

        % welch

        win = hann(W/2);
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

mp_per=mean(beta_per,1);
mp_welch=mean(beta_welch,1);
mp_Yule=mean(beta_Yule,1);
mp_burg=mean(beta_burg,1);
xl=90;


tau_x=1:(xl/dt-W);
tau1=corr(mp_per(1:(xl/dt-W))',tau_x','type','kendall');tau2=corr(mp_welch(1:(xl/dt-W))',tau_x','type','kendall');
tau3=corr(mp_Yule(1:(xl/dt-W))',tau_x','type','kendall');tau4=corr(mp_burg(1:(xl/dt-W))',tau_x','type','kendall');
figure
set(gcf,'position',[10 1 1050 900])
subplot(211)
ax1 = gca;
ax1.Position(2)=ax1.Position(2)-0.05;
set(gca,'Position',ax1.Position)
set(gca,'fontweight','bold','linewidth',3,FontSize=19)
set(gca,'XAxisLocation','top','XTicklabels',[1:2/10:3],'YTick',[]);
xlimits = get(ax1,'XLim');
xinc = (xlimits(2)-xlimits(1))/10;
set(ax1,'XTick',[xlimits(1):xinc:xlimits(2)])
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','left');

for n=1:100
    plot(xt,V(n,:))
    hold on
end
hold on
plot(xt,mean(V),'k','linewidth',3)
xlabel(ax1,'Maximum grazing rate - c','FontWeight','bold',FontSize=24),xlabel(ax2,' '),set(ax2 ,'xticklabels', [])
xline(xl,'--k','linewidth',3)
ylabel('Vegetation biomass - V(t)','FontWeight','bold')
text(-0,13,'(b)','FontWeight','bold',FontSize=24),ylim([0 10.5]),xlim([0 100])
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(ax2,'XTick',[0:10:100])
axis normal


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
text(1,min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)]),'- sliding window -',FontSize=14),text(-0,max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] ),'(c)','FontWeight','bold',FontSize=24)
ylim([min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)])-0.5 max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] )+0.5])
ylabel('Spectral exponent','FontWeight','bold') 
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
legend(['Periodogram ' newline     '\tau = ',num2str(tau1)], ...
    ['Welch          '  newline     '\tau = ',num2str(tau2)], ...
    ['Yule-Walker     ' newline     '\tau = ',num2str(tau3)], ...
    ['Burg           '  newline     '\tau = ',num2str(tau4)] , ...
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

text(0,1.05,'(d)','FontWeight','bold',FontSize=24)
