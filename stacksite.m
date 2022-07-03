% Stack residual spectra for each station after the first t* inversion
% The stacked spectra will be subtracted for the second t* inversion
clear;
tic;

span=21; % span for smoothing
ifplotall=1;ifplotind=0;ifoutput=1;
pors='S';
minnfl=10; % Min number of site effect files to stack
paradir='/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo';
stnlst='/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/input/station.lst';
freqlst='/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/input/allfreq.lst';
% indir=[paradir '/tstar_fcp0.5-20MPa0.27res/resspec'];
indir='~/Work/resspec';
outdir=[paradir '/sitespec'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end
if strcmpi(pors,'P')
    maxfreq=10;
elseif strcmpi(pors,'S')
    maxfreq=4;
end

% Read station list
tmp=importdata(stnlst,' ',1);
stnnm=tmp.textdata(2:end);
stnloc=tmp.data(:,1:3);  % Lon, Lat, Ele
stnnm=sort(stnnm);

% Read all frequencies
allfreq=load(freqlst,'-ascii');

if ifplotall
    figure(3);clf;h2=gcf;
    nrow=8;ncol=8;
end

% Stack for each station
isub=0;
for istn=1:length(stnnm)
% for istn=5
    stn=stnnm{istn};
%     if ~strcmpi(stn,'A01')
%         continue
%     end
    infllst=dir(sprintf('%s/*%sresspec_%s*',indir,upper(pors),stn));nfl=length(infllst);
    if nfl<1
        continue
    end
    isub=isub+1;
    outfl=sprintf('%s/%ssite_%s.spec',outdir,upper(pors),stn);
    
    if ifplotind
        figure(1);clf;h1=gcf;
        subplot(121);hold on;box on;title(sprintf('Site effect: %s %d traces',stn,nfl));
        xlabel('Frequency (Hz)');ylabel('Residual specturm (nm/s)');
        xlim([0 maxfreq]);ylim([-3 3]);
        ax1=gca;set(ax1,'fontsize',24,'Xcolor','k','Ycolor','k');
        subplot(122);hold on;box on;title(['Site effect (%): ' sprintf('%s %d traces',stn,nfl)]);
        xlabel('Frequency (Hz)');ylabel('Residual specturm (%)');
        xlim([0 maxfreq]);ylim([-10 10]);
        ax2=gca;set(ax2,'fontsize',24,'Xcolor','k','Ycolor','k');
    end
    if ifplotall
        subplot(nrow,ncol,isub);hold on;box on;title(sprintf('%s: %d',stn,nfl));
        xlabel('Frequency (Hz)');ylabel('Site eff (nm/s)');
        xlim([0 maxfreq]);ylim([-3 3]);
        ax3=gca;set(ax3,'fontsize',10,'Xcolor','k','Ycolor','k');
    end
   
    allspec=nan(length(allfreq),length(infllst));
    allratio=nan(length(allfreq),length(infllst));
    for ifl=1:nfl
%     for ifl=2:4
        fl=infllst(ifl);
        if strcmp(fl.name(1),'.')
            continue
        end
        flnm=[fl.folder '/' fl.name];
        data=load(flnm,'-ascii');
%         indfreq=find(data(:,1)==allfreq);
        [~,indfreq]=ismember(data(:,1),allfreq);
        allspec(indfreq,ifl)=data(:,2);
        allratio(indfreq,ifl)=data(:,3);
        if ifplotind
            plot(ax1,data(:,1),data(:,2),'color',[0.7 0.7 0.7],'linewidth',1);
            plot(ax2,data(:,1),data(:,3),'color',[0.7 0.7 0.7],'linewidth',1);
        end
        if ifplotall
            plot(ax3,data(:,1),data(:,2),'color',[0.7 0.7 0.7],'linewidth',1);
        end
    end
    avespec=mean(allspec,2,'omitnan');
    errspec=std(allspec,0,2,'omitnan');
    averatio=mean(allratio,2,'omitnan');
    errratio=std(allratio,0,2,'omitnan');
    
    if nfl<minnfl
        smavespec=zeros(length(allfreq),1);smerrspec=smavespec;
        smaveratio=zeros(length(allfreq),1);smerrratio=smavespec;
    else
        smavespec=nan(length(allfreq),1);smerrspec=smavespec;
        smaveratio=nan(length(allfreq),1);smerrratio=smavespec;        
        realind=find(isnan(avespec)==0);
        tmp1=smooth(avespec(realind),span);tmp2=smooth(errspec(realind),span);
        smavespec(realind)=tmp1;smerrspec(realind)=tmp2;
        realind=find(isnan(averatio)==0);
        tmp1=smooth(averatio(realind),span);tmp2=smooth(errratio(realind),span);
        smaveratio(realind)=tmp1;smerrratio(realind)=tmp2;
    end
    outdata=[allfreq smavespec smerrspec smaveratio smerrratio];
    
    if ifoutput
%         dlmwrite(outfl,outdata,'precision',4);
        fid=fopen(outfl,'w');
        fprintf(fid,'%10.4f  %15.8e  %15.8e  %6.2f  %6.2f\n',outdata');
        fclose(fid);
    end
    
    if ifplotind
        plot(ax1,allfreq,smavespec,'b','linewidth',3);
%         plot(ax1,allfreq,savespec,'g','linewidth',2);
        plot(ax1,allfreq,smavespec-smerrspec,'b','linewidth',1);
        plot(ax1,allfreq,smavespec+smerrspec,'b','linewidth',1);
        plot(ax2,allfreq,smaveratio,'b','linewidth',3);
        plot(ax2,allfreq,smaveratio-smerrratio,'b','linewidth',1);
        plot(ax2,allfreq,smaveratio+smerrratio,'b','linewidth',1);
        plot(ax1,allfreq,avespec,'r','linewidth',3);
        plot(ax1,allfreq,avespec-errspec,'r','linewidth',1);
        plot(ax1,allfreq,avespec+errspec,'r','linewidth',1);
        plot(ax2,allfreq,averatio,'r','linewidth',3);
        plot(ax2,allfreq,averatio-errratio,'r','linewidth',1);
        plot(ax2,allfreq,averatio+errratio,'r','linewidth',1);
        set(h1,'renderer','painters','PaperPositionMode','auto');
        print(h1,sprintf('%s/%s_%sresidual.eps',outdir,stn,upper(pors)),'-depsc');
    end
    if ifplotall
        plot(ax3,allfreq,avespec,'r','linewidth',3);
        plot(ax3,allfreq,avespec-errspec,'r','linewidth',1);
        plot(ax3,allfreq,avespec+errspec,'r','linewidth',1);
        plot(ax3,allfreq,smavespec,'b','linewidth',3);
        plot(ax3,allfreq,smavespec-smerrspec,'b','linewidth',1);
        plot(ax3,allfreq,smavespec+smerrspec,'b','linewidth',1);
    end
%     savespec=smavespec;
end
if ifplotall
%     set(h2,'renderer','painters','PaperPositionMode','auto');
%     print(h2,sprintf('%s/all_%sresidual.eps',outdir,upper(pors)),'-depsc');
    print(h2,sprintf('%s/all_%sresidual.png',outdir,upper(pors)),'-dpng');
end
tend=toc;
fprintf('%d min and %.1f s\n',floor(tend/60),rem(tend,60));
