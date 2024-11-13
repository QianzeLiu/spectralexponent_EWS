clear
clc
SLP=ncread('H:\slp.mon.mean.nc','slp');
lon=ncread('H:\slp.mon.mean.nc','lon');
lat=ncread('H:\slp.mon.mean.nc','lat');
time=ncread('H:\slp.mon.mean.nc','time');
t0=datetime(1800,1,1);
T=t0+double(time)/24;
x=length(60:96);y=(length(8:26));
N=831;
S=SLP(60:96,8:26,:);W1=zeros(x*y,N);W2=zeros();


for n=1:x*y
    for k=1:N
        W1(n,k)=S(n+(k-1)*x*y);
    end
end

for n=1:x*y
    for k=0:828/12-1
        W2(n,1+5*k:5+5*k)=W1(n,11+k*12:15+k*12);
    end
end
slpmonths=mean(W2);

xt=1947:1988; slpyears=zeros(1,length(xt));
for k=0:828/12-1
    slpyears(k+1)=mean(mean(W2(:,1+k*5:5+k*5)));
end

XT=1946:(76-46)/(155-1):1976;

step=119;W=step+1;
nfft=W;Fs=W;          
c=0.01;d=0.1;
a=W*c/2;b=W*d/2;
P1=zeros();P2=zeros();P3=zeros();P4=zeros();
for m=1:155-step                      
    Y=slpmonths(m:m+step);



    %Burg
    order = 10;
    [pxx, ~] = pburg(Y, order);
    x =c:(d-c)/(length(a:b)-1):d;
    y =pxx(ceil(a:b));
    y=transpose(y);
    logx = log10(x);
    logy = log10(y);
    beta_burg(m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

    acf=autocorr(Y,1);
    AC1(m)=acf(2);
    VAR(m)=var(Y);

end


%%
figure
set(gcf,'position',[100 100 850 800])
subplot(411)
plot(xt,slpyears(1:length(xt)),'k','linewidth',3),xlim([1946 1988]),ylim([1008,1016])
line([1946 1976],[mean(slpyears(1:length(47:76))),mean(slpyears(1:length(47:76)))],'Color','k','LineStyle','--','linewidth',3)
line([1976 1988],[mean(slpyears(length(47:76)+1:length(xt))) mean(slpyears(length(47:76)+1:length(xt)))],'Color','k','LineStyle','--','linewidth',3)
ylabel('Presure（hPa）','FontWeight','bold')
text(1946,1018,'(a)','FontWeight','bold',FontSize=24)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(gca,'xticklabels', [])

subplot(412)
plot(XT(W:end),beta_burg,'g','linewidth',3),xlim([1946 1988])
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
title('Spectral exponent','FontWeight','bold')
text(1946,-11.2,'(b)','FontWeight','bold',FontSize=24)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(gca,'xticklabels', [])
ylim([min(beta_burg)-.1 max( beta_burg )+.1])


subplot(413)
plot(XT(W:end),AC1,'color','#f2811d','linewidth',3),xlim([1946 1988])
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
title('AC 1','FontWeight','bold')

text(1946,4.1e-4,'(c)','FontWeight','bold',FontSize=24)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
set(gca,'xticklabels', [])
ylim([min(AC1)-.1 max( AC1 )+.1])

subplot(414)
plot(XT(W:end),VAR,'color','#8e6fad','linewidth',3),xlim([1946 1988])
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
xlabel('Year','FontWeight','bold')
title('Variance','FontWeight','bold')
text(1946,6.5,'(d)','FontWeight','bold',FontSize=24)
text(1947,7.7e-3,'────     sliding window     ────',FontSize=14)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
ylim([min(VAR)-0.5 max( VAR )+.5])
