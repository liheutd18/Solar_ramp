clear 
load ace15mins
load solar5
load solar15
x1 = solar_15min;
x2 = solar5min;    % Contains NaN data
x3=reshape(x2, 3, size(x2, 1)/3)';x4=(repmat(x1, 1, 3));
tmp_error_nl=x3-x4;
% Solar_error_max = max(tmp_error_nl, [], 2);
% Solar_error_min = min(tmp_error_nl, [], 2);
x2=0;
y = ACEmin;
%x1=Solar_error_max;x2=Solar_error_min;
%Compute the regression coefficients for a linear model with an interaction term.
X = [ones(size(x1)) x1 x2 x1.*x2];

X1=X(1:2880,:);y1=y(1:2880,:);
b = regress(y1,X1);    % Removes NaN data

YFIT = b(1) + b(2)*X(2881:end,2) + b(3)*X(2881:end,3) + b(4)*X(2881:end,4);
Yreal=y(2881:end);
plot(Yreal)
hold on
plot(YFIT)
%Plot the data and the model.
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):100:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('error_max')
ylabel('error_min')
zlabel('ACE')
view(50,10)
hold off