clear
close all

[cumulative1_INP,outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP,datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP,topmodelsx_INP, M_INP,flag1_INP,typedecline2_INP]=options

datevecfirst1=datevecfirst1_INP;

DT=DT_INP;

cadtemporal='daily';
caddisease='coronavirus'
datatype='cases'
cadregion='USA'
outbreakx=52
caddate1='05-11-2020'

npatches_fixed=4;
onset_fixed=0;
typedecline2=[1 2 3];
smoothfactor1=1;

flag1=1;
method1=0;
dist1=0;
calibrationperiod1=90;
forecastingperiod=30;
topmodels1=1:3;


getperformance=1;

[RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=bagging_forecasts(npatches_fixed,onset_fixed,typedecline2,smoothfactor1,cadtemporal,caddisease,datatype,cadregion,outbreakx,caddate1,flag1,method1,dist1,calibrationperiod1,forecastingperiod,topmodels1,getperformance,DT,datevecfirst1)

figure(101)
subplot(2,2,1)

line1=plot(MSEFS_model1(:,1),MSEFS_model1(:,2),'ko-')
set(line1,'LineWidth',2)
xlabel('Time')
ylabel('MSE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,2)

line1=plot(MAEFS_model1(:,1),MAEFS_model1(:,2),'ko-')
set(line1,'LineWidth',2)
xlabel('Time')
ylabel('MAE')

set(gca,'FontSize', 16);
set(gcf,'color','white')


subplot(2,2,3)

line1=plot(PIFS_model1(:,1),PIFS_model1(:,2),'ko-')
set(line1,'LineWidth',2)
xlabel('Time')
ylabel('Coverage 95% PI')

set(gca,'FontSize', 16);
set(gcf,'color','white')


subplot(2,2,4)

line1=plot(WISFS(:,1),WISFS(:,2),'ko-')
set(line1,'LineWidth',2)
xlabel('Time')
ylabel('WIS')

set(gca,'FontSize', 16);
set(gcf,'color','white')

