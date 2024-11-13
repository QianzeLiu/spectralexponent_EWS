
for n=1:N1
    for m=1:2000-step                       
        Y=Z1(n,m:m+step);

        acf=autocorr(Y,1); 
        AC1(n,m)=acf(2);
        VAR(n,m)=var(Y);
    end
end
     meanAC1=mean(AC1);
 meanVAR=mean(VAR);


%%
xl=100;
xt=10*dt:10*dt:t;
tau_x=1:(xl/dt-W);
tau1=corr(mp_burg(1:(xl/dt-W))',tau_x','type','kendall');tau2=corr(meanAC1(1:(xl/dt-W))',tau_x','type','kendall');
tau3=corr(meanVAR(1:(xl/dt-W))',tau_x','type','kendall');
figure
set(gcf,'position',[10 10 800 1000])
subplot(411)
for n=1:N1
    plot(xt,Z1(n,:))
    hold on
end
xlim([00 700]),ylabel('Z(t)','FontWeight','bold'),text(000,0.9,'(a)','FontWeight','bold',FontSize=24)

hold off
set(gca,'fontweight','bold','linewidth',4,FontSize=19)

subplot(412)
err=std(beta_burg);


plot(xt(W+1:xl/dt),mp_burg(1:(xl/dt-W)),'g','linewidth',2.5)
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[mp3(1:(xl/dt-W))-abs(err(1:(xl/dt-W))),fliplr(mp3(1:(xl/dt-W))+abs(err(1:(xl/dt-W)) ))],'r');
pa.FaceColor = "g";
pa.FaceAlpha = 0.2;  
pa.LineStyle = 'none';  
title('Spectral exponent','FontWeight','bold')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlim([000 700]),
set(gca ,'xticklabels', [])
text(000,2,'(b)','FontWeight','bold',FontSize=24)
text(600,mp_burg(xl/dt-W),['\tau = ',num2str(tau1)],'color','k','FontWeight','bold',FontSize=14)
xline(xt(W+1),'--','color',[0.5 0.5 0.5],'linewidth',3)

subplot(413)
err=std(AC1);


plot(xt(W+1:xl/dt),meanAC1(1:(xl/dt-W)),'color','#f2811d','linewidth',2.5)
title('AC1','FontWeight','bold')
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[meanAC1(1:(xl/dt-W))-abs(err(1:(xl/dt-W))),fliplr(meanAC1(1:(xl/dt-W))+abs(err(1:(xl/dt-W)) ))],'r');
pa.FaceColor = '#f2811d'; 
pa.FaceAlpha = 0.2; 
pa.LineStyle = 'none';  
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlim([0 700])
set(gca ,'xticklabels', [])
text(00,1.3,'(c)','FontWeight','bold',FontSize=24)
text(600,meanAC1(xl/dt-W),['\tau = ',num2str(tau2)],'color','k','FontWeight','bold',FontSize=14)
xline(xt(W+1),'--','color',[0.5 0.5 0.5],'linewidth',3)


subplot(414)
err=std(VAR);


plot(xt(W+1:xl/dt),meanVAR(1:(xl/dt-W)),'color','#8e6fad','linewidth',2.5)
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[meanVAR(1:(xl/dt-W))-abs(err(1:(xl/dt-W))),fliplr(meanVAR(1:(xl/dt-W))+abs(err(1:(xl/dt-W)) ))],'r');
pa.FaceColor = '#8e6fad'; 
pa.FaceAlpha = 0.2;  
pa.LineStyle = 'none';  
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlim([0 700]),
xlabel('Time','FontWeight','bold'),text(150,1.7,'(b)','FontWeight','bold',FontSize=24)
title('Variance','FontWeight','bold')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
text(00,0.005,'(d)','FontWeight','bold',FontSize=24)
axis normal
text(600,meanVAR(xl/dt-W),['\tau = ',num2str(tau3)],'color','k','FontWeight','bold',FontSize=14)
xline(xt(W+1),'--','color',[0.5 0.5 0.5],'linewidth',3)

text(1,meanVAR(1),'-sliding window-',FontSize=12)


%%
for n=1:N1
    tau_SE(n)=corr(beta_burg(n,1:(xl/dt-W))',tau_x','type','kendall');
    tau_AC1(n)=corr(AC1(n,1:(xl/dt-W))',tau_x','type','kendall');
    tau_VAR(n)=corr(VAR(n,1:(xl/dt-W))',tau_x','type','kendall');

end

tau=abs([tau_SE;tau_AC1;tau_VAR;]);
threshold=0.01:0.01:1;
color={'g','#f2811d','#8e6fad'};
num=tau;
figure
set(gcf,'position',[10 100 700 700])
for n=1:3
    for i=1:length(threshold)
        th=threshold(i);
        num(n,i)=length(find(tau(n,:)>=th));
    end

    plot(threshold,num(n,:)/100,'color',cell2mat(color(n)),'linewidth',3)
    hold on
end
auc=sum(0.01.*num/100,2);
legend(['Spectral exponent  AUC = ',num2str(auc(1))], ...
    ['AC 1                         AUC = ',num2str(auc(2))], ...
    ['Variance                  AUC = ',num2str(auc(3))],FontSize=12,Location="southwest")
legend('boxoff')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)

set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlabel('Kendall \tau threshold','FontWeight','bold')
ylabel('Proportion of regime shift detected','FontWeight','bold')
axis tight
text(-0,1.05,'(e)','FontWeight','bold',FontSize=24)