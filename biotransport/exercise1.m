function res =exercise1(initial, time)
if nargin < 2 
    initial = [19/3; 0; 108/100;0;0;0];
%     initial = [19; 0; 110; 0; 0; 0];
end
if nargin < 1
    time = 100;
end
% tspan = linspace(0,time);
[t,vrbls] = ode45(@ODEs, [0 time], initial);
hold on
subplot (3,1,1);
plot (t, 3*vrbls(:,1));
ylabel('Insulin(plasma)', 'FontSize', 20);
title('Constant dose 216mg/min', 'FontSize', 20);
subplot(3,1,2);
plot (t, 11*vrbls(:,2));
ylabel('Insulin(remote)', 'FontSize', 20);
subplot(3,1,3);
plot(t, vrbls(:,3));
xlabel('Time(min)', 'FontSize', 20);
ylabel('Glucose', 'FontSize', 20);


end
function res = ODEs(t, inputs)
I=108;
t1 = 6;
t2=100;
t3=36;
E=0.2;
V1=3;
V2=11;
x = inputs(1);
y = inputs (2);
z = inputs (3);
h1= inputs(4);
h2 = inputs(5);
h3 = inputs(6);
dxdt= f1(z) - E*(x/V1 - y/V2)-x/t1;
dydt = E*(x/V1 - y/V2)-y/t2;
dzdt = f5(h3) +I-f2(z) - f3(z)*f4(y);
dh1dt = 3*(x-h1)/t3;
dh2dt = 3*(h1-h2)/t3;
dh3dt = 3*(h2-h3)/t3;
t
res = [dxdt; dydt; dzdt; dh1dt; dh2dt; dh3dt];

end

function res = f1(z)
if nargin <1
    z = linspace(1,400);
end
V3= 10;
exponent = -100*z./(300*V3)+6.6;
for i =1:length(exponent)
    func1(i) = 209/(1+exp(exponent(i)));
end
res = func1;
end
function res = f2(z)
if nargin <1
    z = linspace(1,400);
end
V3= 10;
exponent = -z.*1/(144/V3);
for i =1:length(exponent)
    func2(i) = 72*(1-exp(exponent(i)));
end
res = func2;
end
function res = f3(z)

if nargin <1
    z = linspace(1,400);
end
V3= 10;
func3 = 1000*.01*z./V3;
set(0, 'defaultaxesfontsize', 18);
 set(0, 'defaulttextfontsize', 18)
 hold on
plot(z, func3);

res = func3;
end
function res = f4(y)

if nargin <1
    y = linspace(1,1000);
end
V2= 11;
t2=100;
E= .2;
nonvariablepart= 1/V2 +1/(E*t2);
for i =1:length(y)
    insidelog = 11*y(i)*nonvariablepart;
    func4(i) = 90/(1+exp(-1.772*log(insidelog)+7.76))+4;
end
res = func4;
end
function res = f5(x)
clf
if nargin <1
    x = linspace(10,200);
end
V1= 3;
for i =1:length(x)
    exponent= 3*0.29*x(i)/V1 -7.5;
    func5(i) =180/(1+exp(exponent)); 
end
res = func5;
end