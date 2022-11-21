function data2=getData(cadtemporal,caddisease,datatype,cadregion,DT,datevecfirst1,datevecend1,date1,outbreak1,forecastingperiod)

filename1=strcat('./input/cumulative-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-',datestr(datenum(datevecend1),'mm-dd-yy'),'.txt');

data=load(filename1);

dataprov=data';

data1=dataprov(outbreak1,:)';

data1=[data1(1);diff(data1)];

datenum1=datenum(datevecfirst1):DT:datenum(date1);


data2=data1(length(datenum1):1:length(datenum1)+forecastingperiod-1);

