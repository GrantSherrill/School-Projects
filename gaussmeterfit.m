clear all

B=[0.000001,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4]';
data=[60,17.42,8.75,6.45,3.58,2.44,1.53,1.26,0.99,0.82,0.74]';


A=linspace(0,10,1000);
umguess=[0.23,-2,1];

modelfun = @(cons) linfun(cons,data,B);
[con] = fminsearch(modelfun,umguess);

bestapprox=con(2)+(con(1)./ (con(3)+(B.^2)));


newline=con(2)+(con(1)./ (con(3)+(A.^2)));

figure()
hold on

% plot(B, bestapprox,'LineWidth',1.8)

plot(B,data,'.','MarkerSize',28)
plot(A, newline,'LineWidth',1.8)
legend('Fit function','Data')
xlabel('Distance(cm)')
ylabel('Field(mT)')
h = gca;
set(h,'FontSize',[24]);
hold off







function J = linfun(c,obs,B)

y=c(2)+c(1)./( c(3)+(B.^2));
J = (y-obs)'*((y-obs));

end 