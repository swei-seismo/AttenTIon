% % % Program to display misfits of all t* inversions with different parameters
% % % Shawn Wei, March 2015
clear;

workdir=('/Users/sowei/Work/Lau/Qtomo/alltstar');
% bandlst=[3];
bandlst=[1 1.5 2]; % fc(P)/fc(S)
alphalst=[0. 0.1 .2 0.27 .3 .4 .5 0.6];
sigmalst=[0.5 5 10 20];
subdblst=[1 2 3 4 11 22 33 44];

i=1;
% for iband=1:length(bandlst)
%     bandwidth=bandlst(iband);
%     banddir=[workdir 'alltstar' num2str(bandwidth)];
    loglst=dir([workdir '/event*_fc*log']);
    for ilog=1:length(loglst)
        logfl=loglst(ilog).name;
        paras=strsplit(logfl,{'_','MPa','sub','.log'});
        stress=char(paras(2));
        if (strfind(stress,'-')>0)
            sigma=0.5;
        else
            sigma=sscanf(stress,'%3s%d');
            sigma=sigma(end);
        end
        if (strfind(stress,'fcp')>0)
            fcps=1;
        elseif (strfind(stress,'fcs')>0)
            fcps=1.5;
        elseif (strfind(stress,'psp')>0)
            fcps=2;
        end
        alpha=str2double(char(paras(3)));
        subdb=str2double(char(paras(4)));
        logfl=[workdir '/' logfl];
%         [~,last2]=system(['tail -n 2 ' logfl]);
%         last2=str2num(last2);
%         if isempty(last2) % sum
%             [~,lastline]=system(['tail -n 1 ' logfl]);
%             lastline=str2num(lastline);
%             resP=lastline(1)/lastline(3);
%             misP=lastline(2)/lastline(3);
%             numP=lastline(3);
%             resS=lastline(4)/lastline(6);
%             misS=lastline(5)/lastline(6);
%             numS=lastline(6);
%         else              % average
            [~,lastline]=system(['tail -n 1 ' logfl]);
            lastline=str2num(lastline);
            resP=lastline(1);
            misP=lastline(2);
            numP=lastline(3);
            resS=lastline(4);
            misS=lastline(5);
            numS=lastline(6);
%         end
        result(i,:)=[fcps alpha sigma subdb resP misP numP resS misS numS];
%         result(i,:)=[fcps alpha sigma subdb misP resP numP misS resS numS];
        i=i+1;
    end
% end
% result=[bandwidth alpha sigma subdb resP misP numP resS misS numS];

figure(3);clf;
bandwd=1;
ind3=find((result(:,1)==bandwd)&(result(:,4)==1));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'rp','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==11));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'rp','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==2));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'r^','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==22));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'r^','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==3));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'ro','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==33));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'ro','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==4));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'rs','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==44));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'rs','MarkerSize',20);hold on;
bandwd=1.5;
ind3=find((result(:,1)==bandwd)&(result(:,4)==1));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'bp','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==11));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'bp','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==2));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'b^','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==22));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'b^','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==3));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'bo','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==33));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'bo','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==4));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'bs','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==44));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'bs','MarkerSize',10);hold on;
bandwd=2;
ind3=find((result(:,1)==bandwd)&(result(:,4)==1));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'kp','MarkerFaceColor','k','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==11));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'kp','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==2));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'k^','MarkerFaceColor','k','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==22));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'k^','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==3));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'ko','MarkerFaceColor','k','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==33));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'ko','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==4));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'ks','MarkerFaceColor','k','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==44));
plot3(result(ind3,2),result(ind3,3),result(ind3,5),'ks','MarkerSize',10);hold on;
box on;grid on;
h=legend('subdb1-1.0','subdb11-1.0','subdb2-1.0','subdb22-1.0','subdb3-1.0','subdb33-1.0',...
    'subdb4-1.0','subdb44-1.0',...
    'subdb1-1.5','subdb11-1.5','subdb2-1.5','subdb22-1.5','subdb3-1.5','subdb33-1.5',...
    'subdb4-1.5','subdb44-1.5');
h.FontSize=16;
xlabel('\alpha');ylabel('Stress Drop (MPa)'),zlabel('Average misfit');
title(['t*(P)']);
set(gca,'FontSize',16);

figure(4);clf;
bandwd=1;
ind3=find((result(:,1)==bandwd)&(result(:,4)==1));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'rp','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==11));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'rp','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==2));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'r^','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==22));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'r^','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==3));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'ro','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==33));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'ro','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==4));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'rs','MarkerFaceColor','r','MarkerSize',20);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==44));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'rs','MarkerSize',20);hold on;
bandwd=1.5;
ind3=find((result(:,1)==bandwd)&(result(:,4)==1));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'bp','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==11));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'bp','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==2));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'b^','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==22));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'b^','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==3));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'bo','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==33));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'bo','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==4));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'bs','MarkerFaceColor','b','MarkerSize',10);hold on;
ind3=find((result(:,1)==bandwd)&(result(:,4)==44));
plot3(result(ind3,2),result(ind3,3),result(ind3,8),'bs','MarkerSize',10);hold on;
% bandwd=2;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==1));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'kp','MarkerFaceColor','k','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==11));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'kp','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==2));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'k^','MarkerFaceColor','k','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==22));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'k^','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==3));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'ko','MarkerFaceColor','k','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==33));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'ko','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==4));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'ks','MarkerFaceColor','k','MarkerSize',10);hold on;
% ind3=find((result(:,1)==bandwd)&(result(:,4)==44));
% plot3(result(ind3,2),result(ind3,3),result(ind3,8),'ks','MarkerSize',10);hold on;
box on;grid on;
h=legend('subdb1-1.0','subdb11-1.0','subdb2-1.0','subdb22-1.0','subdb3-1.0','subdb33-1.0',...
    'subdb4-1.0','subdb44-1.0',...
    'subdb1-1.5','subdb11-1.5','subdb2-1.5','subdb22-1.5','subdb3-1.5','subdb33-1.5',...
    'subdb4-1.5','subdb44-1.5');
h.FontSize=16;
xlabel('\alpha');ylabel('Stress Drop (MPa)'),zlabel('Average misfit');
title(['t*(S)']);
set(gca,'FontSize',16);

% All datasets
% % result(i,:)=[fcps alpha sigma subdb resP misP numP resS misS numS];
allresult=zeros(length(bandlst)*length(alphalst)*length(sigmalst),size(result,2)-1);
j=1;
for iband=1:length(bandlst)
    for ialpha=1:length(alphalst)
        for isigma=1:length(sigmalst)
            for i=1:size(result,1) % Sum up all datasets for same parameters
%                if ((result(i,1)==bandlst(iband)) && (result(i,4)>10) && ...
%                        (result(i,2)==alphalst(ialpha)) && (result(i,3)==sigmalst(isigma)))
               if ((result(i,1)==bandlst(iband)) && ...
                       (result(i,2)==alphalst(ialpha)) && (result(i,3)==sigmalst(isigma)))
                   allresult(j,4)=allresult(j,4)+result(i,5)*result(i,7);
                   allresult(j,5)=allresult(j,5)+result(i,6)*result(i,7);
                   allresult(j,6)=allresult(j,6)+result(i,7);
                   allresult(j,7)=allresult(j,7)+result(i,8)*result(i,10);
                   allresult(j,8)=allresult(j,8)+result(i,9)*result(i,10);
                   allresult(j,9)=allresult(j,9)+result(i,10);
               end
            end
            allresult(j,1:3)=[bandlst(iband) alphalst(ialpha) sigmalst(isigma)];
            allresult(j,4)=allresult(j,4)/allresult(j,6);
            allresult(j,5)=allresult(j,5)/allresult(j,6);
            allresult(j,7)=allresult(j,7)/allresult(j,9);
            allresult(j,8)=allresult(j,8)/allresult(j,9);
            j=j+1;
        end
    end
end
% % % allresult=[bandwidth alpha sigma resP misP numP resS misS numS];

figure(1);clf;
ind3=find(allresult(:,1)==1);
plot3(allresult(ind3,2),allresult(ind3,3),allresult(ind3,4),'r*','MarkerSize',10);hold on;
ind4=find(allresult(:,1)==1.5);
plot3(allresult(ind4,2),allresult(ind4,3),allresult(ind4,4),'b^','MarkerSize',10);hold on;
% ind4=find(allresult(:,1)==2);
% plot3(allresult(ind4,2),allresult(ind4,3),allresult(ind4,4),'k^','MarkerSize',10,'MarkerFaceColor','k');hold on;
box on;grid on;
h=legend('fc(S)=fc(P)','fc(S)=fc(P)/1.5','dt*(S-P)');
h.FontSize=16;
xlabel('\alpha');ylabel('Stress Drop (MPa)');zlabel('Average misfit');
title('t*(P)');
set(gca,'FontSize',16);

figure(2);clf;
ind3=find(allresult(:,1)==1);
plot3(allresult(ind3,2),allresult(ind3,3),allresult(ind3,7),'r*','MarkerSize',10);hold on;
ind4=find(allresult(:,1)==1.5);
plot3(allresult(ind4,2),allresult(ind4,3),allresult(ind4,7),'b^','MarkerSize',10);hold on;
% ind4=find(allresult(:,1)==2);
% plot3(allresult(ind4,2),allresult(ind4,3),allresult(ind4,7),'k^','MarkerSize',10,'MarkerFaceColor','k');hold on;
box on;grid on;
h=legend('fc(S)=fc(P)','fc(S)=fc(P)/1.5','dt*(S-P)');
h.FontSize=16;
xlabel('\alpha');ylabel('Stress Drop (MPa)');zlabel('Average misfit');
title('t*(S)');
set(gca,'FontSize',16);

figure(5);clf;
ind3=find(allresult(:,1)==1 & allresult(:,3)==0.5);
plot(allresult(ind3,2),allresult(ind3,4),'rp','MarkerSize',20,'MarkerFaceColor','r');hold on;
ind3=find(allresult(:,1)==1 & allresult(:,3)==5);
plot(allresult(ind3,2),allresult(ind3,4),'b^','MarkerSize',20,'MarkerFaceColor','b');hold on;
ind3=find(allresult(:,1)==1 & allresult(:,3)==10);
plot(allresult(ind3,2),allresult(ind3,4),'ro','MarkerSize',20,'MarkerFaceColor','r');hold on;
ind3=find(allresult(:,1)==1 & allresult(:,3)==20);
plot(allresult(ind3,2),allresult(ind3,4),'bs','MarkerSize',20,'MarkerFaceColor','b');hold on;
box on;grid on;
h=legend('\Delta\sigma = 0.5-20 MPa','\Delta\sigma = 5 MPa',...
    '\Delta\sigma = 10 MPa','\Delta\sigma = 20 MPa');
h.FontSize=24;
xlabel('\alpha');ylabel('Average misfit');xlim([-0.1 0.7]);
title('t*(P)');
set(gca,'FontSize',24);
