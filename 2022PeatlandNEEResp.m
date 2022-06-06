%%%%%%%% analyse monthly peatland NEE responses to temperature %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (manuel.helbig@dal.ca) on 2022-06-06
% used for Helbig et al. (2022) Warming response of peatland CO2 sink is
% sensitive to seasonality in warming trends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
path = '/Users/manuelhelbig/MEGA/GWF_Postdoc/SPRUCE/';

%% read list of sites
data = importdata('Table_SITES.csv');
P_MAP=data.data(:,6); % annual precipitation
TA_MAP=data.data(:,4); % annual air temperature
lat=data.data(:,1); % latitude
lon=data.data(:,2); % longitude

[P_MAP, i] = sort(P_MAP,'ascend');
TA_MAP=TA_MAP(i);
lat=lat(i);
lon=lon(i);

alphabet    = 'a':'t'
[map] = brewermap(20,'PuBu');
for k=1:20
    if lon(k) <-20
    plot(TA_MAP(k),P_MAP(k),'s','MarkerFaceColor',map(k,:));
    hold on
    axis square
    text(TA_MAP(k)+0.2,P_MAP(k)-20,alphabet(k))

    elseif lon(k) > -20 & lon(k) < 35
    plot(TA_MAP(k),P_MAP(k),'o','MarkerFaceColor',map(k,:));
    text(TA_MAP(k)+0.2,P_MAP(k)-20,alphabet(k))
    hold on
    axis square
    else

    plot(TA_MAP(k),P_MAP(k),'^','MarkerFaceColor',map(k,:));
    hold on
    axis square
    text(TA_MAP(k)+0.2,P_MAP(k)-20,alphabet(k))
    end
end

%% read monthly net ecosystem exchange data
data=importdata([path '2022PeatlandNEE_Dataset.csv']);

NEE=data.data(:,6);     % net ecosystem exchange
MNT=data.data(:,5);     % month
YY=data.data(:,4);      % year
SITE=data.data(:,1);    % site ID
LAT=data.data(:,2);     % latitude
LON=data.data(:,3);     % longitude
TA=data.data(:,7);      % air temperature
TA_CRU=data.data(:,11); % air temperature (Climate Research Unit)
WTD=data.data(:,8);     % water table depth
SWIN=data.data(:,9);    % shortwave incoming radiation
EVI_MOD=data.data(:,10);% Enhanced Vegetation Index (MODIS)
PRE=data.data(:,12);    % precipitation (Climate Research Unit)
PET=data.data(:,13);    % potential evapotranspiration (Climate Research Unit)
TS=data.data(:,14);     % near-surface soil temperature
TS_DP=data.data(:,15);  % deep soil temperature
VPD=data.data(:,16);    % vapour pressure deficit (local)
VPD_CRU=data.data(:,17);% vapour pressure deficit (Climate Research Unit)
SD=data.data(:,18);     % snow depth (local)

LABEL=strsplit(data.textdata{1},','); % variable names
uSITE=unique(SITE);     % site IDs

%% classification of seasons
for n=1:24 % loop through sites
    mn=MNT(SITE==n);
    ta=TA_CRU(SITE==n);
    lat(n)=median(LAT(SITE==n));
    for l=1:12
        % get monthly averages
        if ~isempty(mn)
            TAmean(n,l)=nanmean(ta(mn==l));
        else
            TAmean(n,l)=NaN;
        end
    end
end
% standardise seasonality
for n=1:24;
    TAmean(n,:)=(TAmean(n,:)-mean(TAmean(n,:)))./std(TAmean(n,:));
end

for n=1:24;
    if ~isnan(sum(TAmean(n,:)))
        SPR(n)=find(TAmean(n,:)>0,1,'first');
        FAL(n)=find(TAmean(n,:)>0,1,'last');
        PEAK(n)=find(TAmean(n,:)==max(TAmean(n,:)));
        MINM(n)=find(TAmean(n,:)==min(TAmean(n,:)));
    end
end
SPR(SPR==0)=NaN;
FAL(FAL==0)=NaN;
PEAK(PEAK==0)=NaN;
MINM(MINM==0)=NaN;

%% calculate standard deviation for variables (per site & month)
for k=1:length(uSITE)
    for n=1:12
        ind=find(SITE==uSITE(k) & MNT==n);
        x=TA(ind);
        y=NEE(ind);
        z=WTD(ind);
        w=TS(ind);
        if sum(~isnan(y))>=3
            stdNEE(k,n)=nanstd(y);
            stdTA(k,n)=nanstd(x);
            stdWTD(k,n)=nanstd(z);
            stdTS(k,n)=nanstd(w);
        else
            stdNEE(k,n)=NaN;
            stdTA(k,n)=NaN;
            stdWTD(k,n)=NaN;
            stdTS(k,n)=NaN;
        end
        clear x y z w
    end
end
%% plot standard deviation (NEE & Ta)
figure,
subplot(1,2,1)
plotshaded(1:12,[nanmean(stdNEE)'-(nanstd(stdNEE)./sqrt(sum(~isnan(stdNEE))))' nanmean(stdNEE)'+(nanstd(stdNEE)./sqrt(sum(~isnan(stdNEE))))'],'r');
hold on
plot(1:12,nanmean(stdNEE));
axis square
xlim([0.5 12.5])
ylabel('\sigma NEE [g C m^{-2} mnt^{-1}]')

subplot(1,2,2)
plotshaded(1:12,[nanmean(stdTA)'-(nanstd(stdTA)./sqrt(sum(~isnan(stdTA))))' nanmean(stdTA)'+(nanstd(stdTA)./sqrt(sum(~isnan(stdTA))))'],'r');
hold on
l1=plot(1:12,nanmean(stdTA));

plotshaded(1:12,[nanmean(stdTS)'-(nanstd(stdTS)./sqrt(sum(~isnan(stdTS))))' nanmean(stdTS)'+(nanstd(stdTS)./sqrt(sum(~isnan(stdTS))))'],'r');
hold on
l2=plot(1:12,nanmean(stdTS));
legend([l1 l2],'air temperature','soil temperature');
axis square
xlim([0.5 12.5])
ylabel('\sigma T [\circC]')

%% calculate anomaly time series 

aNEE=NaN(size(NEE));
aTA=NaN(size(NEE));
aTS=NaN(size(NEE));
aTA_CRU=NaN(size(NEE));
aSW=NaN(size(NEE));
aWTD=NaN(size(WTD));
aEVI_MOD=NaN(size(NEE));
aPRE=NaN(size(NEE));

for s=1:12
    ind=find(MNT==s);
    SIT_MNT=SITE(ind);
    uSIT_MNT=unique(SIT_MNT);
    for r=1:length(uSIT_MNT)
        % get anomalies of variables per month and site
        i=find(MNT==s & SITE==uSIT_MNT(r));
        aNEE(i)=(NEE(i)-nanmean(NEE(i)));
        aTA_CRU(i)=(TA_CRU(i)-nanmean(TA_CRU(i))); 
        aTA(i)=(TA(i)-nanmean(TA(i)));
        aTS(i)=TS(i)-nanmean(TS(i));
        aWTD(i)=WTD(i)-nanmean(WTD(i));
        aSW(i)=SWIN(i)-nanmean(SWIN(i));
        aPRE(i)=PRE(i)-nanmean(PRE(i));
        aEVI_MOD(i)=EVI_MOD(i)-nanmean(EVI_MOD(i));
    end
end

%% linear mixed effects models for each month
clear SD SE
for s=1:12
    
    ind=find(MNT==s);
    x = aTA(ind); % code can be applied to other variables (change variable here)
    y = aNEE(ind);
    z = SITE(ind);
    n_length(s)=sum(~isnan(y) & ~isnan(x));
    tbl = table(x,y,z,'VariableNames',{'Ta','NEE','Site'});  
    lme = fitlme(tbl,'NEE~Ta+(1|Site)+(Ta-1|Site)');
    [beta(:,s),betanames,stats] = fixedEffects(lme);
    s
    pval(s)=double(stats(2,6));
    CI(s,:)=stats(2,7:8);

    SD(s)=double(stats(2,3)).*sqrt(double(stats(2,5)));
    SE(s)=double(stats(2,3));
end

%% plot monthly regression slopes
CI=double(CI);
neg=beta(2,:)-CI(:,1)';
pos=beta(2,:)-CI(:,2)';

figure,
errorbar(1:12,beta(2,:),neg,pos);
hold on

refline(0,0)
xlim([0.5 12.5]);
clear CI beta neg pos
% add significance levels
for k=1:12
    text(k-0.4,3.7,['\it{n} = ' num2str(n_length(k))]);
    if pval(k)<0.001
        text(k-0.5,3.6,['***']);
    elseif pval(k)<0.01
        
        text(k-0.5,3.6,['**']);
    elseif  pval(k)<0.05
        text(k-0.5,3.5,['*']);
    end
end
ylabel('Sensitivity_{T_a - NEE} [g C m^{-2} mnt^{-1} \circC^{-1}]');

%% add season definition
SEASONS=NaN(size(aNEE));
for n=1:24
    SEASONS(MNT==SPR(n) & SITE == n) = 2; % spring
    SEASONS(MNT==FAL(n) & SITE == n) = 5; % fall
    SEASONS(MNT<FAL(n) & SITE == n & MNT>=PEAK(n)) = 4; % late summer
    SEASONS(MNT>SPR(n) & SITE == n & MNT<PEAK(n)) = 3; % early summer
    if MINM(n)<10
       SEASONS(MNT<SPR(n) & SITE == n & MNT>=MINM(n)) = 1; % late winter
    else
       SEASONS((MNT<SPR(n) | MNT>=MINM(n)) & SITE == n) = 1; % late winter  
    end
    if MINM(n)<10
       SEASONS((MNT>FAL(n) | MNT<MINM(n)) & SITE == n) = 6; % early winter 
    else
       SEASONS(MNT>FAL(n) & MNT<MINM(n) & SITE == n) = 6; % early winter  
    end
end

%% derive NEE sensitivity to TA for each season
clear R2
clear beta2 pval2 NEEbin CI2 pos2 neg2 TEMPbin2 n_length2 NEEbin2 n_length TEMPbin
for s=1:length(unique(SEASONS))
    ind=find(SEASONS==s);
    x = aTA(ind); % replace variable here to run for other variables
    y = aNEE(ind);
    z = SITE(ind);
    
    tbl = table(x,y,z,'VariableNames',{'Ta','NEE','Site'});
    lme = fitlme(tbl,'NEE~Ta+(1|Site)+(Ta-1|Site)');
    [beta2(:,s),betanames,stats] = fixedEffects(lme);
    
    NEEbin(s)=nanmean(NEE(ind)); % mean NEE for season
    n_length(s)=sum(~isnan(y) & ~isnan(x));
    pval2(s)=double(stats(2,6));
    CI2(s,:)=stats(2,7:8);
    clear ind
end

CI2=double(CI2);
neg1=beta2(2,:)-CI2(:,1)';
pos1=beta2(2,:)-CI2(:,2)';
neg_T=neg1;
pos_T=pos1;

figure, % plot sensitivities

SEAS=categorical({'late winter';'spring';'early summer';'late summer';'fall';'early winter'});
SEAS = reordercats(SEAS,{'late winter';'spring';'early summer';'late summer';'fall';'early winter'});

errorbar(SEAS,beta2(2,:),neg1,pos1);
hold on
plot(SEAS(pval2<0.05),beta2(2,pval2<0.05),'o');
scatter(SEAS,beta2(2,:),25,NEEbin,'filled');
caxis([-30 30]);
colormap(flipud(brewermap(8,'BrBg')))
ylabel('Sensitivity_{T_a-NEE} [g C m^{-2} mnt^{-1} \circC^{-1}]');

xlim=get(gca,'xlim');
plot(xlim,[0 0],'k');

h_bar = colorbar('location','NorthOutside','FontSize',14);
set(get(h_bar,'ylabel'),'String', {'NEE [g C m^{-2} mnt^{-1}]'},'FontSize',14);
axis square
xtickangle(45)
for r=1:length(unique(SEASONS));
    ht=text(SEAS(r),0.28,num2str(n_length(r)));
    set(ht,'Rotation',45)
end