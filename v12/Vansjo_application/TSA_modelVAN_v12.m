% Script to run MyLake (v_12) for Vansjo 
% Used in E-project Pcode
% by TSA, last modified 23.02.2005

clear all;
path(path,'H:\Mylake_v12_package\air_sea') %path for air-sea toolbox
path(path,'H:\Mylake_v12_package\v12') %path for MyLake model code

global ies80 Eevapor;
test_time=0;
Eevapor=0;

load H:\Mylake_v12_package\Vansjo_application\Observations\vansjoice.dat
load H:\Mylake_v12_package\Vansjo_application\Observations\vansjotemp.dat 

[Obs_TP_Chla, trash]=xlsread('H:\Mylake_v12_package\Vansjo_application\Observations\Vansjo_TP_Chla84_05.xls');
% year, month, day, TP, Chla
load H:\Mylake_v12_package\Vansjo_application\Observations\Obs_SS.dat; 
% year, month, day, SS (mg/L)
load H:\Mylake_v12_package\Vansjo_application\Observations\Obs_PO4P.dat; 
% year, month, day, PO4P (microg/L)

lake='Vansjo-Storefjorden';
year=1994;
m_start=[1994,5,1]; 
m_stop=[2000,12,31];
%m_stop=[2005,11,30];

initfile='H:\Mylake_v12_package\Vansjo_application\VAN_init_v12.xls';
parafile='H:\Mylake_v12_package\Vansjo_application\VAN_para_v12.xls';
inputfile='H:\Mylake_v12_package\Vansjo_application\VAN_input_PGS84_00_v12.xls';

tic
        [zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
        P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt]...
           = solvemodel_v12(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake');    
run_time=toc

DoF_realtime=DoF+tt(1)-1; %TSA, antatt at tidsteg er 1 dag
DoM_realtime=DoM+tt(1)-1; %TSA
DoF_plottime=DoF+tt(1)-1-datenum(year,1,1); %TSA, antatt at tidsteg er 1 dag
DoM_plottime=DoM+tt(1)-1-datenum(year,1,1); %TSA

tt_mod = tt - datenum(year,1,1); %time now scaled so that it begins from the 1 january of the "year" (=0)

%=Ice thickness observations (tt_mod, Hice (m))
IceObs=[datenum(fliplr(vansjoice(:,1:3))) - datenum(year,1,1),vansjoice(:,4)./100];
inx=find((IceObs(:,1)<(datenum(m_start)- datenum(year,1,1)))|(IceObs(:,1)>(datenum(m_stop)- datenum(year,1,1))));
IceObs(inx,:)=[];

%=Temperature profile observations (tt_mod, z, T)
Datestr=num2str(vansjotemp(:,1));
Dummydate=[str2num(Datestr(:,1:4)),str2num(Datestr(:,5:6)),str2num(Datestr(:,7:8))];
TempObs=[datenum(Dummydate) - datenum(year,1,1), -vansjotemp(:,3)./100, vansjotemp(:,4)];

%=align temperature observations with model results
alku=[1;find(diff(TempObs(:,1))~=0)+1];
loppu=[find(diff(TempObs(:,1))~=0); length(TempObs)];

for i=1:length(alku)
inxt=find(tt_mod==TempObs(alku(i),1));
    if (isempty(inxt)==0)
    TempMod(alku(i):loppu(i))=interp1(zz,Tzt(:,inxt),TempObs(alku(i):loppu(i),2));
    else
    TempMod(alku(i):loppu(i))=NaN;    
    end    
end

zlim = [0 max(zz)];
tlim = [min(tt_mod) max(tt_mod)];

% thermocline depth
zt = MixStat(12,:);
 

figure(1)
clf
subplot(312)
pcolor(tt_mod,zz,Tzt)
shading interp
axis ij
hold on
plot(tt_mod,zt,'k','LineWidth',1);
hold off
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0 24]);
colorbar;
set(gca,'fontsize',9);
H=get(gca,'Position');
ylabel('Depth (m)')
set(gca,'TickDir','out')
%set(gca,'xticklabel',xlabs)
%colormap(colmap)

subplot(311);
plot(tt_mod,His(1,:)+His(2,:),'-r',tt_mod,His(1,:)-His(3,:),':c',tt_mod,His(1,:),'-b')
hold on
plot(IceObs(:,1),IceObs(:,2),'m*')
H(2)=H(2)+0.23; %0.226
set(gca,'Position',H);
set(gca,'fontsize',9);
set(gca,'ylim',[0 1]);
datetick('x','mmm');
set(gca,'XTickLabel',[]);
set(gca,'YTick',[0.2 0.4 0.6 0.8 1]);
set(gca,'TickDir','out')
ylabel('Hice, Hsnow (m)');

grid on;


figure(2)
clf
plot([0 25],[0 25],':b', TempObs(:,3), TempMod, '.r');
set(gca,'fontsize',9);
ylabel('Modelled temperature (^oC)');
xlabel('Measured temperature (^oC)');
axis([0 25 0 25]);
axis square;
title([lake ' ' datestr(datenum(m_start),28) '--' datestr(datenum(m_stop),28)]); 
grid on;

figure(3)
clf
fign=6;
inxpos=find(isnan(TempMod(alku))==0);
posalku=alku(inxpos);
posloppu=loppu(inxpos);

for i = 1:min(fign,length(posalku))
   subplot(3,3,i);
   inxt=find(tt_mod==TempObs(posalku(i),1));
   plot(Tzt(:,inxt),zz,'-b',TempObs(posalku(i):posloppu(i),3),TempObs(posalku(i):posloppu(i),2),'.r');
   axis([0 25 zlim])
   axis ij
   title(['Vansjø-Stfj. ' datestr(tt_mod(inxt)+datenum(year,1,1)),],'fontsize',8); 
   set(gca,'FontSize',8) 
end;

   subplot(334)
   ylabel('Depth (m)','fontsize',8)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
   subplot(331)
   ylabel('Depth (m)','fontsize',8)
    
    subplot(335)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
     subplot(336)
   xlabel('Temperature (^{o}C)','fontsize',8)

figure(4)
clf
for i = fign+1:min(12,length(posalku))
   subplot(3,3,i-fign);
   set(gca,'FontSize',8) 
   inxt=find(tt_mod==TempObs(posalku(i),1));
   plot(Tzt(:,inxt),zz,'-b',TempObs(posalku(i):posloppu(i),3),TempObs(posalku(i):posloppu(i),2),'.r');
   axis([0 25 zlim])
   axis ij
   title(['Vansjø-Stfj. ' datestr(tt_mod(inxt)+datenum(year,1,1)),],'fontsize',8); 
end;

   subplot(334)
   ylabel('Depth (m)','fontsize',8)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
    subplot(331)
   ylabel('Depth (m)','fontsize',8)
   
    subplot(335)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
     subplot(336)
   xlabel('Temperature (^{o}C)','fontsize',8)


figure(5)
clf
subplot(2,2,1);
pcolor(tt_mod,zz,Czt)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 1]);
colorbar;
set(gca,'fontsize',9);
title(['Dissolved passive tracer']); 


subplot(2,2,2);
pcolor(tt_mod,zz,Szt)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0 0.01]);
colorbar;
set(gca,'fontsize',9);
title(['Sedimenting passive tracer (g/L)']); 

 subplot(223)
 plot(Az./1e+6,zz)
 axis ij
 xlabel('Horizontal area (km^2)')
 ylabel('Depth (m)')
 set(gca,'Ylim',[0 40])
 grid on

 subplot(224)
 plot(tt_mod,MixStat(4,:),'g')
 hold on
 plot(tt_mod,-MixStat(5,:),'r')
 datetick('x','mmm');
 ylabel('Mean growth/loss rates (1/day)')
 grid on
 title([lake ' ' datestr(datenum(m_start),28) '--' datestr(datenum(m_stop),28)]); 

 
figure(6)
clf
subplot(3,1,1);
pcolor(tt_mod,zz,Pzt)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
%plot(tt_mod,MixStat(3,:),'w-')
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0 30]);
colorbar;
set(gca,'fontsize',9);
title(['Dissolved phosphorus (\mug/l)']); 

subplot(3,1,2);
pcolor(tt_mod,zz,Chlzt)
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
%plot(tt_mod,MixStat(3,:),'w-')
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0 15]);
colorbar;
set(gca,'fontsize',9);
title(['Chlorophyll {\ita} (\mug/l)']); 

subplot(3,1,3);
pcolor(tt_mod,zz,PPzt)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
%plot(tt_mod,MixStat(3,:),'w-')
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0 15]);
colorbar;
set(gca,'fontsize',9);
title(['Particulate inorg. phosphorus (\mug/l)']); 


figure(7)
clf
subplot(321);
pcolor(tt_mod,zz,squeeze(P3zt_sed(:,:,1)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
%caxis([0 1]);
colorbar;
set(gca,'fontsize',8);
title('Diss. P in sediment porewater')
ylabel('mg m^{-3}'); 

subplot(323);
pcolor(tt_mod,zz,squeeze(P3zt_sed(:,:,2)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
%caxis([0 1]);
colorbar;
set(gca,'fontsize',8);
title('P in inorg. sediment solids')
ylabel('mg kg^{-1}'); 

subplot(325);
pcolor(tt_mod,zz,(1/2.5)*(1-squeeze(P3zt_sed(:,:,4))).*squeeze(P3zt_sed(:,:,3)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
%caxis([0 1]);
colorbar;
set(gca,'fontsize',8);
title('Chl {\ita} in sediment solids')
ylabel('mg kg^{-1}'); 

subplot(322);
pcolor(tt_mod,zz,squeeze(P3zt_sed_sc(:,:,1)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
%caxis([0  0.1]);
colorbar;
set(gca,'fontsize',8);
title('Diss. P source from sediment')
ylabel('mg/m^{-3} day^{-1}'); 

subplot(324);
pcolor(tt_mod,zz,squeeze(P3zt_sed_sc(:,:,2)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0  0.5]);
colorbar;
set(gca,'fontsize',8);
title('Inorg. partic. P source from sediment')
ylabel('mg/m^{-3} day^{-1}'); 

subplot(326);
pcolor(tt_mod,zz,squeeze(P3zt_sed_sc(:,:,3)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([0 0.5]);
colorbar;
set(gca,'fontsize',8);
title('Chl {\ita} source from sediment')
ylabel('mg/m^{-3} day^{-1}'); 


figure(8)
clf
zinx=find(zz<4); %depth layer considered
kk=0;
dum=datevec(tt_mod+datenum(year,1,1));
yrs=dum(:,1);
months=dum(:,2);

monthly_qual=NaN*ones(12*(max(yrs)-min(yrs)+1),12);
F_OM=1e+6*0.012;    %mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)
for i=min(yrs):max(yrs)
   for j=1:12
    kk=kk+1;   
    tinx=find((yrs==i)&(months==j)); 

    monthly_qual(kk,1)=i; %year
    monthly_qual(kk,2)=j; %month
    monthly_qual(kk,3)=15; %day
    monthly_qual(kk,4)=datenum(i,j,15); %dtenumber
    monthly_qual(kk,5)=mean(mean(Pzt(zinx,tinx)+PPzt(zinx,tinx)+DOPzt(zinx,tinx)+Chlzt(zinx,tinx))); %TotP
    monthly_qual(kk,6)=mean(mean(Chlzt(zinx,tinx))); %Chla
    monthly_qual(kk,7)=mean(mean(Szt(zinx,tinx)+Chlzt(zinx,tinx)./F_OM)); %SS
    monthly_qual(kk,8)=mean(mean(Chlzt(zinx,tinx))); %Chla
    monthly_qual(kk,9)=mean(mean(Pzt(zinx,tinx))); %P
    monthly_qual(kk,10)=mean(mean(PPzt(zinx,tinx))); %PP
    monthly_qual(kk,11)=mean(mean(DOPzt(zinx,tinx))); %DOP
       zdum=-9.6 + 3.6*log10(Pzt(zinx,tinx)+PPzt(zinx,tinx)+DOPzt(zinx,tinx)+Chlzt(zinx,tinx)) + 0.23*Tzt(zinx,tinx);
       if(isempty(zdum)==0)
        monthly_qual(kk,12)=max(mean(1./(1+exp(-zdum)))); %max(P[>10% cyanobacteria])
       end
   end
 end
 
 %synchronize observations to model simulations
         TP_Chla_Mod=NaN*ones(length(Obs_TP_Chla),2); 
        for i=1:length(Obs_TP_Chla)
            inx=find(datenum(Obs_TP_Chla(i,1:3))==monthly_qual(:,4));
            if(isempty(inx)==0)
                TP_Chla_Mod(i,1)=monthly_qual(inx,5); %1) TP
                TP_Chla_Mod(i,2)=monthly_qual(inx,6); %2) Chla
            end
        end
 
 yearly_qual=NaN*ones((max(yrs)-min(yrs)+1),6);
 kk=0;
for i=min(yrs):max(yrs)
    kk=kk+1;   
    tinx=find((monthly_qual(:,1)==i)&(monthly_qual(:,2)>5)&(monthly_qual(:,2)<10)); 
    yearly_qual(kk,1)=i; %year
    yearly_qual(kk,2)=mean(monthly_qual(tinx,5)); %TotP
    yearly_qual(kk,3)=mean(monthly_qual(tinx,6)); %Chla
    yearly_qual(kk,4)=mean(monthly_qual(tinx,7)); %SS
    yearly_qual(kk,5)=mean(monthly_qual(tinx,10)); %PIP
    yearly_qual(kk,6)=max(monthly_qual(tinx,12)); %P[cyano>10%]
 end
 
subplot(311)
 plot(monthly_qual(:,4),monthly_qual(:,5),'.-')
 hold on
 plot(datenum(Obs_TP_Chla(:,1:3)),Obs_TP_Chla(:,4),'r.')
 set(gca,'ylim',[0 70]);
 ylabel('TotP  (\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 
 subplot(312)
 plot(monthly_qual(:,4),monthly_qual(:,6),'.-')
 hold on
 plot(datenum(Obs_TP_Chla(:,1:3)),Obs_TP_Chla(:,5),'r.')
 set(gca,'ylim',[0 20]);
 ylabel('Chl a  (\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 
 subplot(313)
 plot(monthly_qual(:,4),monthly_qual(:,9),'.-')
 hold on
 plot(datenum(Obs_PO4P(:,1:3)),Obs_PO4P(:,4),'r.') %mg/m3
 set(gca,'ylim',[0 40]);
 ylabel('PO4 (\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 
figure(9)
clf
for i=1:length(Obs_TP_Chla)
inx=find(datenum(Obs_TP_Chla(i,1:3))==monthly_qual(:,4));
 if(isempty(inx)==0)
 subplot(221)  
 plot(Obs_TP_Chla(i,4),monthly_qual(inx,5),'b.')
 hold on
 subplot(222)  
 plot(Obs_TP_Chla(i,5),monthly_qual(inx,6),'b.')
 hold on
 end
end
subplot(221)
axis('square')
axis([0 60 0 60])
xlabel('Observed')
ylabel('Modelled')
plot([0 60],[0 60],'--')
title(['TotP ({\mu}g/l)']);

subplot(222)
axis('square')
axis([0 20 0 20])
xlabel('Observed')
ylabel('Modelled')
plot([0 30],[0 30],'--')
title(['Chla {\ita} ({\mu}g/l)']); 

subplot(223)
for i=1:length(Obs_SS)
inx=find(datenum(Obs_SS(i,1:3))==monthly_qual(:,4));
 if(isempty(inx)==0)
 loglog(1e-3*Obs_SS(i,4),monthly_qual(inx,7),'b.') %mg/L -> g/L
 hold on
 end
end
axis('square')
axis([1e-3 1e-1 1e-3 1e-1])
xlabel('Observed')
ylabel('Modelled')
loglog([1e-5 1],[1e-5 1],'--')
title(['Susp. matter {\ita} (g/l)']); 

subplot(224)
for i=1:length(Obs_PO4P)
inx=find(datenum(Obs_PO4P(i,1:3))==monthly_qual(:,4));
 if(isempty(inx)==0)
 plot(Obs_PO4P(i,4),monthly_qual(inx,9),'b.') %mg/m3
 hold on
 end
end
axis('square')
axis([0 20 0 20])
xlabel('Observed')
ylabel('Modelled')
plot([0 30],[0 30],'--')
title(['PO4 (\mug/l)']); 

figure(10)
clf
subplot(221)
plot(tt_mod,MixStat(3,:))
 ylabel('Total PAR extinction coeff. at 2 m  (1/m)')
 set(gca,'ylim',[1 2]);
 grid on
 datetick('x','mmm');
 title([lake ' ' datestr(datenum(m_start),28) '--' datestr(datenum(m_stop),28)]); 

subplot(222);
pcolor(tt_mod,zz,Qzt_sed)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
datetick('x','mmm');
set(gca,'ylim',zlim);
%axis([tlim zlim]);
caxis([-3 3]);
colorbar;
set(gca,'fontsize',9);
title(['Heat flux from sediment (W/m^2)']); 

subplot(223);
plot((P3zt_sed(:,end,5)-P3zt_sed(:,end,6))./Az,zz)
title('Net sedimentation during simulation')
xlabel('kg/m^2')
axis ij


subplot(224);
 plot(P3zt_sed(:,end,5)-P3zt_sed(:,end,6),zz)
title('Net sedimentation during simulation')
xlabel('kg')
axis ij
 

figure(11)
clf
subplot(3,1,1);
pcolor(tt_mod,zz,Szt)
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 0.05]);
colorbar;
set(gca,'fontsize',9);
title(['Susp. matter (g/l)']); 

subplot(3,1,2);
pcolor(tt_mod,zz,PPzt)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 15]);
colorbar;
set(gca,'fontsize',9);
title(['Particulate inorg. phosphorus (\mug/l)']); 

subplot(3,1,3);
pcolor(tt_mod,zz,DOPzt)
axis ij
shading interp
hold on
plot(DoF_plottime,zeros(1,length(DoF_plottime)),'wv','MarkerSize',5);
plot(DoM_plottime,zeros(1,length(DoM_plottime)),'w^','MarkerSize',5);
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 15]);
colorbar;
set(gca,'fontsize',9);
title(['Dissolved org. phosphorus (\mug/l)']); 

figure(12)
clf

subplot(221);
pcolor(tt_mod,zz,DOCzt)
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 3000]);
colorbar;
set(gca,'fontsize',9);
title('DOC (mg/m^{-3})'); 


subplot(222);
pcolor(tt_mod,zz,DOPzt)
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 50]);
colorbar;
set(gca,'fontsize',9);
title('DOP (mg/m^{-3})'); 

subplot(223);
pcolor(tt_mod,zz,squeeze(P3zt_sed(:,:,4)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0.8 1]);
colorbar;
set(gca,'fontsize',8);
title('Fraction inorg. matter in sediment solids')
ylabel('m^3 m^{-3}'); 

subplot(224);
pcolor(tt_mod,zz,(Chlzt./F_OM)./(Szt+Chlzt./F_OM))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 0.3]);
colorbar;
set(gca,'fontsize',8);
title('Fraction organic matter in susp solids')
ylabel('kg^3 kg^{-3}'); 

figure(13)
clf
subplot(221);
pcolor(tt_mod,zz,1e-3*squeeze(P3zt_sed(:,:,5)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
colorbar;
set(gca,'fontsize',8);
title('Sedimentation')
ylabel('tons'); 

subplot(222);
pcolor(tt_mod,zz,1e-3*squeeze(P3zt_sed(:,:,6)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
colorbar;
set(gca,'fontsize',8);
title('Resuspension')
ylabel('tons'); 

subplot(223);
pcolor(tt_mod,zz,squeeze(P3zt_sed(:,:,7)))
shading interp
axis ij
datetick('x','mmm');
set(gca,'ylim',zlim);
caxis([0 1e-6/0.05]);
colorbar;
set(gca,'fontsize',8);
title('Fraction of new sediment')
 
figure(14)
clf
subplot(211)
 plot(monthly_qual(:,4),monthly_qual(:,8),'r-')
 hold on
 plot(monthly_qual(:,4),sum(monthly_qual(:,8:9)'),'b-')
 plot(monthly_qual(:,4),sum(monthly_qual(:,8:10)'),'m-')
 plot(monthly_qual(:,4),sum(monthly_qual(:,8:11)'),'k-')
 set(gca,'ylim',[0 60]);
 ylabel('(\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 H=legend('Chla','dissP','PIP','DOP');
 set(H,'fontsize',8);
 
 subplot(212)
 plot(monthly_qual(:,4),monthly_qual(:,7),'.-')
 hold on
 plot(datenum(Obs_SS(:,1:3)),1e-3*Obs_SS(:,4),'r.:') %mg/L -> g/L
 set(gca,'ylim',[0 0.032]);
 ylabel('Susp. matter  (g/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 
 figure(15)
 clf
 subplot(211)
 plot(tt_mod, cumsum(MixStat(13,:)),'r')
 hold on
 plot(tt_mod, cumsum(MixStat(14,:)),'b')
 plot(tt_mod, cumsum(MixStat(15,:)),'c')
 plot(tt_mod, (MixStat(18,:)),'m')
 plot(tt_mod, cumsum(MixStat(13,:)-MixStat(14,:)-MixStat(15,:)+MixStat(16,:)+MixStat(17,:))-MixStat(18,:),'k')
 
 datetick('x','mmm');
 H=legend('Inflow of P', 'Outflow of P', 'Sedimentation of P', 'Change in lake P', 'P-Balance');
 set(H,'fontsize',8);
 grid on;
  ylabel('kg')
  
 dum=datevec(tt_mod+datenum(year,1,1));
 yrs=dum(:,1);
 kk=0;
 Intern=NaN*ones(length(max(yrs)-min(yrs))+1,5);
 for i=min(yrs):max(yrs)
 inx=find(yrs==i);
 kk=kk+1;
 Intern(kk,1)=sum(MixStat(16,inx));
 Intern(kk,2)=sum(MixStat(17,inx));
 Intern(kk,3)=sum(MixStat(19,inx)); %Net flux out of sediment
 Intern(kk,4)=sum(MixStat(13,inx)); %Inflow
 Intern(kk,5)=sum(MixStat(14,inx)); %Outflow
 Intern(kk,6)=sum(MixStat(20,inx)); %Algae available inflow
 end    
 
 subplot(212)
 plot([min(yrs)+1:max(yrs)], -1e-3*Intern(2:end,3),'k','linewidth',2)
 hold on
 plot([min(yrs)+1:max(yrs)], 1e-3*Intern(2:end,4),'r','linewidth',2)
 plot([min(yrs)+1:max(yrs)], 1e-3*Intern(2:end,6),'r--','linewidth',2)
 plot([min(yrs)+1:max(yrs)], 1e-3*Intern(2:end,5),'b','linewidth',2)

 H=legend('Net TP flux to sed.','TP inflow','P04 inflow','TP outflow');
 title('P budget','fontweight','bold')
 set(H,'fontsize',8);
 ylabel('tons/year')
 grid on;
 
 
 %=Figures for modelpaper
 tlims=[datenum([1984,12,15]) datenum([2001,1,15])];
 
 figure(20)
 clf
 subplot(311)
 plot(monthly_qual(:,4),monthly_qual(:,5),'.-')
 hold on
 plot(datenum(Obs_TP_Chla(:,1:3)),Obs_TP_Chla(:,4),'r+')
 set(gca,'ylim',[0 70]);
 ylabel('mg/m^3','fontsize',9)
 title('TotP (monthly mean, 0-4m)','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
 set(gca,'xlim',tlims)
 
 subplot(312)
 plot(monthly_qual(:,4),monthly_qual(:,6),'.-')
 hold on
 plot(datenum(Obs_TP_Chla(:,1:3)),Obs_TP_Chla(:,5),'r+')
 set(gca,'ylim',[0 20]);
 ylabel('mg/m^3','fontsize',9)
 title('P_{Chla} (monthly mean, 0-4m)','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
 set(gca,'xlim',tlims)
 
 subplot(313)
 %zdum=-9.6 + 3.6*log10(Pzt(zinx,:)+PPzt(zinx,:)+DOPzt(zinx,:)+Chlzt(zinx,:)) + 0.23*Tzt(zinx,:);
 %plot(tt_mod+datenum(year,1,1),mean(1./(1+exp(-zdum))),'r--'); %P[>10% cyanobacteria] daily mean
 %hold on
 bar(yearly_qual(2:end,1),yearly_qual(2:end,6),'c') %yearly mean
 set(gca,'ylim',[0 1]);
 ylabel('Probability','fontsize',9)
 title('P[>10% cyanobacteria] (yearly maximum, 0-4m)','fontweight','bold')
 %datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 %set(gca,'xticklabel',[]);
 set(gca,'xlim',[1984 2001])
 xlabel('year')
 set(gca,'tickdir','out')
 
 figure(21)
 clf
 subplot(311)
 plot(monthly_qual(:,4),monthly_qual(:,9),'.-')
 hold on
 plot(datenum(Obs_PO4P(:,1:3)),Obs_PO4P(:,4),'r+') %mg/m3
 set(gca,'ylim',[0 40]);
 ylabel('mg/m^3','fontsize',9)
 title('P_D (monthly mean, 0-4m)','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
 set(gca,'xlim',tlims)
 
 subplot(312)
 plot(monthly_qual(:,4),monthly_qual(:,7),'.-')
 hold on
 plot(datenum(Obs_SS(:,1:3)),1e-3*Obs_SS(:,4),'r+') %mg/L -> g/L
 set(gca,'ylim',[0 0.032]);
 ylabel('kg/m^3','fontsize',9)
 xlabel('year','fontsize',9)
 title('Susp. matter (monthly mean, 0-4m)','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xlim',tlims)

 
figure(22)
clf
subplot(311)
inx=find(round(TempObs(:,2))<2);
H=plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt(1,:),'-');
 set(gca,'ylim',[0 25]);
 ylabel('^oC','fontsize',9)
 title('Temperature  0-1m','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
set(gca,'xlim',tlims)
 
subplot(312)
inx=find((round(TempObs(:,2))==10)|(round(TempObs(:,2))==11));
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt(11,:),'-');
 set(gca,'ylim',[0 25]);
 ylabel('^oC','fontsize',9)
 title('Temperature  10-11m','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
 set(gca,'xlim',tlims)
 
 subplot(313)
inx=find((round(TempObs(:,2))==30)|(round(TempObs(:,2))==31));
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt(31,:),'-');
 set(gca,'ylim',[0 25]);
 ylabel('^oC','fontsize',9)
 xlabel('year','fontsize',9)
 title('Temperature 30-31m','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xlim',tlims)
 %===============
 
 disp(['Sum of P sinks: '  num2str(round(sum(MixStat(14,:)+MixStat(15,:)))) ' kg']); 
 disp(['Sum of P sources: '  num2str(round( sum(sum(Intern)) + sum(MixStat(13,:)) )) ' kg']);
 
 disp(['Average summer season TotP, Chla, and SS: ' num2str(mean(yearly_qual(:,2:4)))])

 

 