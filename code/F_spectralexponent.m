clc
clear
dt=0.5;t=1000;N=t/dt;
N1=100;xt=dt:dt:t;
Z=zeros(N1,N);
for n=1:N1
    for i=2:N
        Z(n,i)=Z(n,i-1)-dt*(4*(Z(n,i-1))^3+2*(3-xt(i-1)/200)*Z(n,i-1))+dt*0.05*randn/sqrt(dt);
    end
end



%%
step=399;W=400;
nfft=W;         
c=0.01;d=0.1;
a=W*c/2;b=W*d/2;fs=2;
xl=70;tau_x=1:(xl/dt-W);

for n=1:N1
    for m=1:2000-step                       
        Y=Z(n,m:m+step);
        %periodogram
        Fs=length(Y);
        window=boxcar(length(Y));
        [pxx,f]=periodogram(Y,window,nfft,fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_per(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

        %welch
        N = length(Y);
        win = hann(N/2); 
        [pxx, f] = pwelch(Y, win, [], [], fs);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_welch(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);


        % Yule-Walker
        order = 10; 
        [pxx, f] = pyulear(Y, order);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_yule(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

        % Burg
        order = 10; 
        [pxx, f] = pburg(Y, order);
        x =c:(d-c)/(length(a:b)-1):d;
        y =pxx(a:b);
        y=transpose(y);
        logx = log10(x);
        logy = log10(y);
        beta_burg(n,m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

    end


end
%%
mp_per=-mean(beta_per);mp_welch=-mean(beta_welch);mp_Yule=-mean(beta_yule);mp_burg=-mean(beta_burg);


xt=10*dt:10*dt:t;


tau1=corr(mp_per(1:(xl/dt-W))',tau_x','type','kendall');tau2=corr(mp_welch(1:(xl/dt-W))',tau_x','type','kendall');
tau3=corr(mp_Yule(1:(xl/dt-W))',tau_x','type','kendall');tau4=corr(mp_burg(1:(xl/dt-W))',tau_x','type','kendall');

figure
set(gcf,'position',[10 1 1050 900])
subplot(211)
for n=1:N1
    plot(xt,Z1(n,:))
    hold on
end
xlim([000 700]),ylabel('Z(t)','FontWeight','bold'),text(000,0.7,'(a)','FontWeight','bold',FontSize=24)
xline(600,'--k','linewidth',3),set(gca,'xticklabels', [])%xticks([])
hold off
ax1 = gca;
ax1.Position(2)=ax1.Position(2)-0.05;
set(gca,'Position',ax1.Position)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
subplot(212)

plot(0,0,'w')
hold on
 err_per=std(beta_per);

 err_welch=std(beta_welch);

 err_Yule=std(beta_yule);

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
pa3.FaceAlpha = 0.5; 
 pa3.LineStyle = 'none';  
hold on
h=plot(xt(W+1:xl/dt),mp_per(1:(xl/dt-W)),'k',xt(W+1:xl/dt),mp_welch(1:(xl/dt-W)),'r',...
    xt(W+1:xl/dt),mp_Yule(1:(xl/dt-W)),'b',xt(W+1:xl/dt),mp_burg(1:(xl/dt-W)),'g','linewidth',3);
xline(600,'--k','linewidth',3)
xlim([000 700])

xlabel('Time','FontWeight','bold')

xline(xt(W+1),'--','color',[0.5 0.5 0.5],'linewidth',3)
text(5,min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)])-0.2,'— — sliding window — —',FontSize=14),
text(-0,max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] )+0.5,'(b)','FontWeight','bold',FontSize=24)
ylim([min([mp_per(1),mp_welch(1),mp_Yule(1),mp_burg(1)])-0.5 max( [mp_per(xl/dt-W),mp_welch(xl/dt-W),mp_Yule(xl/dt-W),mp_burg(xl/dt-W)] )+0.5])
ylabel('Spectral exponent','FontWeight','bold') 
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
legend(h,['Periodogram ' newline  '\tau = ',num2str(tau1)], ...
    ['Welch           ' newline  '\tau = ',num2str(tau2)]', ...
    ['Yule-Walker     ' newline  '\tau = ',num2str(tau3)], ...
    ['Burg             ' newline  '\tau = ',num2str(tau4)] , ...
    FontSize=12,Location='northwest');
legend('boxoff')
axis normal



%%


for n=1:N1
    tau_per(n)=corr(beta_per(n,1:(xl/dt-W))',tau_x','type','kendall');tau_wel(n)=corr(beta_welch(n,1:(xl/dt-W))',tau_x','type','kendall');
    tau_yul(n)=corr(beta_yule(n,1:(xl/dt-W))',tau_x','type','kendall');tau_bur(n)=corr(beta_burg(n,1:(xl/dt-W))',tau_x','type','kendall');
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
text(-0,1.05,'(c)','FontWeight','bold',FontSize=24)
