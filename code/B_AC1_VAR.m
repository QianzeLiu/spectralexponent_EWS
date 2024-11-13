
for n=1:N1
    for m=1:N-step
        Y=B(n,m:m+step);

        acf=autocorr(Y,1);
        AC1(n,m)=acf(2);
        VAR(n,m)=var(Y);
    end
end
meanAC1=mean(AC1);
meanVAR=mean(VAR);


%%
xl=100;
tau_x=1:(xl/dt-W);
tau1=corr(mp_burg(1:(xl/dt-W))',tau_x','type','kendall');tau2=corr(meanAC1(1:(xl/dt-W))',tau_x','type','kendall');
tau3=corr(meanVAR(1:(xl/dt-W))',tau_x','type','kendall');
figure
set(gcf,'position',[10 10 800 1000])
subplot(411)

ax1 = gca;
set(gca,'Position',ax1.Position)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(gca,'XAxisLocation','top','YAxisLocation','right','XDir','reverse','xlim',[0 10],'XTicklabels',[0.1  0.7  1.16  1.5   1.6  1.7 ],'YTick',[]);
xlimits = get(ax1,'XLim');
xinc = (xlimits(2)-xlimits(1))/5;
set(ax1,'XTick',[xlimits(1):xinc:xlimits(2)])

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','left');
for n=1:100
    plot(xt,B(n,:))
    hold on
end
hold on
plot(xt,mean(B),'k','linewidth',3),xlabel(ax1,'Rainfall rate - R','FontWeight','bold'),xlabel(ax2,' '),set(ax2 ,'xticklabels', [])
xline(xl,'--k','linewidth',3)
ylabel('B(t)','FontWeight','bold') ,xlim([0 100]),ylim([-0.5 7.9]),
text(0,12,'(a)','FontWeight','bold',FontSize=24)
set(ax2,'XTick',[0:20:100])
axis normal
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
subplot(412)
err_burg=std(beta_burg);

plot(xt(W+1:xl/dt),mp_burg(1:(xl/dt-W)),'g','linewidth',2.5)
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[mp_burg(1:(xl/dt-W))-abs(err_burg(1:(xl/dt-W))),fliplr(mp_burg(1:(xl/dt-W))+abs(err_burg(1:(xl/dt-W)) ))],'r');
pa.FaceColor = "g";
pa.FaceAlpha = 0.2; 
pa.LineStyle = 'none'; 
title('Spectral exponent','FontWeight','bold')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlim([0 100]),
xline(xl,'--k','linewidth',3);xline(xt(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
set(gca ,'xticklabels', [])
text(-0,1.4,'(b)','FontWeight','bold',FontSize=24),xlim([0 100]),
text(70,mp_burg(xl/dt-W),['\tau = ',num2str(tau1)],'color','k','FontWeight','bold',FontSize=14)


subplot(413)
err_AC1=std(AC1);

plot(xt(W+1:xl/dt),meanAC1(1:(xl/dt-W)),'color','#f2811d','linewidth',2.5)
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[meanAC1(1:(xl/dt-W))-abs(err_AC1(1:(xl/dt-W))),fliplr(meanAC1(1:(xl/dt-W))+abs(err_AC1(1:(xl/dt-W)) ))],'r');
pa.FaceColor = '#f2811d'; 
pa.FaceAlpha = 0.2;  
pa.LineStyle = 'none';  
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlim([0 100])
xline(xl,'--k','linewidth',3);xline(xt(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
set(gca ,'xticklabels', [])
title('AC1','FontWeight','bold')
text(-0,1.04,'(c)','FontWeight','bold',FontSize=24),xlim([0 100]),
text(70,meanAC1(xl/dt-W),['\tau = ',num2str(tau2)],'color','k','FontWeight','bold',FontSize=14)


subplot(414)
err_VAR=std(VAR);
plot(xt(W+1:xl/dt),meanVAR(1:(xl/dt-W)),'color','#8e6fad','linewidth',2.5)
pa = patch([xt(W+1:xl/dt),fliplr(xt(W+1:xl/dt))],[meanVAR(1:(xl/dt-W))-abs(err_VAR(1:(xl/dt-W))),fliplr(meanVAR(1:(xl/dt-W))+abs(err_VAR(1:(xl/dt-W)) ))],'r');
pa.FaceColor = '#8e6fad';
pa.FaceAlpha = 0.2;  
pa.LineStyle = 'none';  
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
xlim([0 100])
xline(xl,'--k','linewidth',3);xline(xt(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
xlabel('Time','FontWeight','bold'),text(150,1.7,'(b)','FontWeight','bold',FontSize=24)
title('Variance','FontWeight','bold')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
text(-0,0.8,'(d)','FontWeight','bold',FontSize=24),xlim([0 100])
text(1,meanVAR(1),'- sliding window -',FontSize=12)
axis normal
text(70,meanVAR(xl/dt-W),['\tau = ',num2str(tau3)],'color','k','FontWeight','bold',FontSize=14)





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