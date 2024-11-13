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

sm=mean(reshape(slpmonths,5,69),2);
slpmonths_anomoly=reshape(reshape(slpmonths,5,69)-sm,1,345);

xt=1947:1988; slpyears=zeros(1,length(xt));
for k=0:828/12-1
    slpyears(k+1)=mean(mean(W2(:,1+k*5:5+k*5)));
end

XT=1946:(76-46)/(155-1):1976;

step=99;W=step+1;
nfft=W;Fs=2;
c=0.01;d=0.1;
a=W*c/2;b=W*d/2;
P1=zeros();P2=zeros();P3=zeros();P4=zeros();
for m=1:155-step
    Y=slpmonths(m:m+step);
    %periodogram
    window=boxcar(length(Y));
    [pxx,f]=periodogram(Y,window,nfft,Fs);
    x =f(find(f>0.01&f<0.1));
    y =pxx(find(f>0.01&f<0.1));
    logx = log10(x);
    logy = log10(y);
    beta_per(m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2+eps);

    % welch
    win = hann(W/2);
    [pxx, f] = pwelch(Y, win, [], [], Fs);
    x =f(find(f>0.01&f<0.1));
    y =pxx(find(f>0.01&f<0.1));
    y=transpose(y);
    logx = log10(x);
    logy = log10(y);
    beta_welch(m)  = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

    %Yule-Walker
    order = 10;
    [pxx, f] = pyulear(Y, order);
    x =f(find(f>0.01&f<0.1));
    y =pxx(find(f>0.01&f<0.1));
    logx = log10(x);
    logy = log10(y);
    beta_Yule(m)  = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

    % Burg
    order = 10;
    [pxx, f] = pburg(Y, order);
    x =f(find(f>0.01&f<0.1));
    y =pxx(find(f>0.01&f<0.1));
    logx = log10(x);
    logy = log10(y);
    beta_burg(m) = sum((logx - mean(logx)).*(logy - mean(logy))) / sum((logx - mean(logx)).^2);

end


figure
set(gcf,'position',[100 100 850 800])
subplot(311)
plot(xt,slpyears(1:length(xt)),'k','linewidth',3),xlim([1946 1988]),ylim([1008,1016])
line([1946 1976],[mean(slpyears(1:length(47:76))),mean(slpyears(1:length(47:76)))],'Color','k','LineStyle','--','linewidth',3)
line([1976 1988],[mean(slpyears(length(47:76)+1:length(xt))) mean(slpyears(length(47:76)+1:length(xt)))],'Color','k','LineStyle','--','linewidth',3)
ylabel('Presure（hPa）','FontWeight','bold')
text(1946,1017,'(a)','FontWeight','bold',FontSize=24)
set(gca ,'xticklabels', [])
set(gca,'fontweight','bold','linewidth',4,FontSize=19)

subplot(312)
plot(XT(W:end),beta_per,'k','linewidth',3),xlim([1946 1988])
% ylim([-3e-4,3e-4])
text(1946,max(beta_per),'(b)','FontWeight','bold',FontSize=24)
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
ylim([min(beta_per)-.2 max( beta_per )+.2])
set(gca ,'xticklabels', [])


subplot(313)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
plot(XT(W:end),beta_burg,'G','linewidth',3),xlim([1946 1988])
text(1946,max(beta_burg),'(c)','FontWeight','bold',FontSize=24)
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
xlabel('Year','FontWeight','bold')
ylabel('Spectral exponent','FontWeight','bold')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
ylim([min(beta_burg)-.1 max( beta_burg )+.1])
text(1947,max(beta_burg),'────     sliding window     ────',FontSize=14)
%%
figure
set(gcf,'position',[100 100 850 800])
subplot(511)
plot(xt,slpyears(1:length(xt)),'k','linewidth',3),xlim([1946 1988]),ylim([1008,1016])
line([1946 1976],[mean(slpyears(1:length(47:76))),mean(slpyears(1:length(47:76)))],'Color','k','LineStyle','--','linewidth',3)
line([1976 1988],[mean(slpyears(length(47:76)+1:length(xt))) mean(slpyears(length(47:76)+1:length(xt)))],'Color','k','LineStyle','--','linewidth',3)
ylabel('Presure（hPa）','FontWeight','bold')
text(1946,1017,'(a)','FontWeight','bold',FontSize=24)
set(gca ,'xticklabels', [])
set(gca,'fontweight','bold','linewidth',4,FontSize=19)

subplot(512)
plot(XT(W:end),beta_per,'k','linewidth',3),xlim([1946 1988])
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
set(gca ,'xticklabels', [])
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
text(1946,max( beta_per ),'(b)','FontWeight','bold',FontSize=24)
title('Periodogram  Method','FontWeight','bold')
ylim([min(beta_per)-.1 max( beta_per )+.1])
subplot(513)
plot(XT(W:end),beta_welch,'r','linewidth',3),xlim([1946 1988])
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
set(gca ,'xticklabels', [])
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
text(1946,max( beta_welch ),'(c)','FontWeight','bold',FontSize=24)
ylabel('Spectral exponent','FontWeight','bold')
title('Welch  Method','FontWeight','bold')
ylim([min(beta_welch)-.0001 max( beta_welch )+.0001])


subplot(514)
plot(XT(W:end),beta_Yule,'b','linewidth',3),xlim([1946 1988])
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
set(gca ,'xticklabels', [])
title('Yule-Walker Method','FontWeight','bold')
ylim([min(beta_Yule)-.001 max( beta_Yule )+.001])
text(1946,max( beta_Yule ),'(d)','FontWeight','bold',FontSize=24)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)


subplot(515)
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
plot(XT(W:end),beta_burg,'G','linewidth',3),xlim([1946 1988])
text(1946,max(beta_burg),'(e)','FontWeight','bold',FontSize=24)
xline(1976,'--k','linewidth',3),xline(XT(W),'--','color',[0.5 0.5 0.5],'linewidth',3)
xlabel('Year','FontWeight','bold')
title('Burg Method','FontWeight','bold')
set(gca,'fontweight','bold','linewidth',4,FontSize=19)
ylim([min(beta_burg)-.1 max( beta_burg )+.1])
text(1947,max(beta_burg),'────     sliding window     ────',FontSize=14)


set(gca,'fontweight','bold','linewidth',4,FontSize=19)
