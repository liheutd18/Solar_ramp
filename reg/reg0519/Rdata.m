%new data
T_reg_4sec1920 = readtable('C:\Users\lxh180005\Desktop\solar ramp\ACE_REG_2019_2020_received 2020_10_31_TRANSFORMED FOR JOSEPHINE.xlsx');
Reg_Feb=T_reg_4sec1920{1:end,7};
ACE_STAR4sec=T_reg_4sec{2:end,2};
REUP4sec=T_reg_4sec{2:end,4};
REUD4sec=T_reg_4sec{2:end,5};
ACE4sec=T_reg_4sec{2:end,3};
k_5_15=225;%75 or 225
a=length(ACE_STAR4sec);a1=reshape(ACE_STAR4sec,k_5_15,a/k_5_15);ACE_STARmin=sum(a1)'/k_5_15;   
a2=reshape(REUP4sec,k_5_15,a/k_5_15);REUPmin=sum(a2)'/k_5_15;   
a3=reshape(REUD4sec,k_5_15,a/k_5_15);REUDmin=sum(a3)'/k_5_15;   
a4=reshape(ACE4sec,k_5_15,a/k_5_15);ACEmin=sum(a4)'/k_5_15;
RegPro = readtable('C:\Users\lxh180005\Desktop\solar ramp\OASIS Data\20200201_20200301_AS_REQ_DAM_RegUp+Dn_20210529_03_10_20_v1_Josephine Transformed.xlsx');

for i=1:720
    RegPro2(4*(i-1)+1,:)=RegPro1(i,:);
    RegPro2(4*(i-1)+2,:)=RegPro1(i,:);
    RegPro2(4*(i-1)+3,:)=RegPro1(i,:);
    RegPro2(4*(i-1)+4,:)=RegPro1(i,:);
end

RegPro3=RegPro2(end-296:end,:);
exceedbaseline=YTest'>RegPro3(:,2);
c1=sum(exceedbaseline);c2=sum(RegPro3(:,2));
a=4;
YPred1=YPred*a;

tmp_error_nl = reshape(YPred1(1:296), 4, length(YPred1(1:296))/4);
error_max = max(tmp_error_nl, [], 1);
error_min = min(tmp_error_nl, [], 1);error_min=-abs(error_min);
Error=[error_max' error_min'];
for i=1:74
    Error1(4*(i-1)+1,:)=Error(i,:);
    Error1(4*(i-1)+2,:)=Error(i,:);
    Error1(4*(i-1)+3,:)=Error(i,:);
    Error1(4*(i-1)+4,:)=Error(i,:);
end
Error1(297,:)=Error1(296,:);
exceednew=YTest'>Error1(:,1);
b1=sum(exceednew);b2=sum(Error1(:,1));
plot(Error1)


t=1:1:120000;
x=ACE_STAR4sec(1:120000);
[up1,lo1] = envelope(x,900,'peak');
param_small = {'Color',[0.9 0.4 0.1],'Linewidth',2};
plot(t,x)
hold on
p1 = plot(t,up1,param_small{:});
plot(t,lo1,param_small{:});
