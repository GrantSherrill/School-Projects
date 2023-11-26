%% Manual Tracking inputs
clear all
filename='Results from MGS_anti2a1_2022-02-22_2.csv';

mu=1.0016*10^(-3)*10^-12; %In N/m^2*s    *10^-12(m^2/um^2)
r=1; %in um

fr=1; %frames per second

data=importdata(filename);
time=data.data(:,2);
bead=data.data(:,1);
xpos=data.data(:,3);
ypos=data.data(:,4);

xpos=xpos/12.3644; %From pixels to microns
ypos=ypos/12.3644;

time=time/fr;
dt=1/fr;

beadinfo=size(bead);
k=1;
j=1;
for k=1:length(bead)
    j=1;
    for i=1:beadinfo(1)
        if bead(i)==k
            beadx(j,k)=xpos(i);
            beady(j,k)=ypos(i);
            beancount(k)=j;
            j=j+1;
            
        end
    end
end 

for k=1:length(beancount)
      tempbeadx=beadx(1:beancount(k)-1,k);
      tempbeady=beady(1:beancount(k)-1,k);
      tempdisp=sqrt(tempbeadx.^2+tempbeady.^2);





  
    for j=1:beancount(k)-1
        newpos(j,k)=(tempbeady(j)-tempbeady(1))^2;
        rms(j,k)=sum(newpos(1:j,k))/j;
    end

    rmstempx=rms(1:beancount(k)-1,k);

    umguess=[1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    modelfun = @(um) linfun2(um , rmstempx,time);
    [um_opt] = fminsearch(modelfun,umguess);
    slope(k)=um_opt(1);

    
    bestapprox(1:length(rmstempx),k) = (time.*um_opt(1));

    umguess2=[1,1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    modelfun2 = @(um2) quadfun2(um2 , rmstempx,time);
    [um_opt2] = fminsearch(modelfun2,umguess2);
    slope2(k)=um_opt2(1);
    diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));
    if um_opt2(2)<0
        diffforce(k)=-imag(diffforce(k));
    end    
    bestapprox2(1:length(rmstempx),k) = (time.*um_opt2(1))+(time.^2 .* um_opt2(2));

end



%% 

for j=1:beancount(k)-1
    newpos(j,k)=(tempbeadx(j)-tempbeadx(1))^2;
    rms(j,k)=sum(newpos(1:j,k))/j;
end

rmstempx=rms(1:beancount(k)-1,k);

umguess=[1];
time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
time=time';
modelfun = @(um) linfun2(um , rmstempx,time);
[um_opt] = fminsearch(modelfun,umguess);
slope(k)=um_opt(1);


bestapprox(1:length(rmstempx),k) = (time.*um_opt(1));

umguess2=[1,1];
time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
time=time';
modelfun2 = @(um2) quadfun2(um2 , rmstempx,time);
[um_opt2] = fminsearch(modelfun2,umguess2);
slope2(k)=um_opt2(1);
diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));
if um_opt2(2)<0
    diffforce(k)=-imag(diffforce(k));
end    
bestapprox2(1:length(rmstempx),k) = (time.*um_opt2(1))+(time.^2 .* um_opt2(2));

%% Bead Starting Position and Force origin videos

BOI=beadx(:,11);
tb=dt*linspace(1,length(BOI),length(BOI));

starts=2;
off=1;
for k=1:starts
    for j=1:length(BOI)-1
        posm(j,k)=(BOI(j)-BOI(k))^2;
        rmsnew(j,k)=sum(posm(1:j,k))/j;
    end
end
rmsnew(1:k,:)=[];

for k=1:starts
    rmsrow=rmsnew(:,k);
    umguess=[1];
    time=dt*linspace(1,length(rmsnew),length(rmsnew));
    time=time';
    modelfun = @(um) linfun2(um , rmsrow,time);
    [um_opt] = fminsearch(modelfun,umguess);
    slope(k)=um_opt(1);
    
    
    bestappr(1:length(rmsrow),k) = (time.*um_opt(1));
    
    umguess2=[1,1];
    time=dt*linspace(1,length(rmsnew),length(rmsnew));
    time=time';
    modelfun2 = @(um2) quadfun2(um2 , rmsrow,time);
    [um_opt2] = fminsearch(modelfun2,umguess2);
    slope2(k)=um_opt2(1);
    diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));
    if um_opt2(2)<0
        diffforce(k)=-imag(diffforce(k));
    end    
    bestappr2(1:length(rmsrow),k) = (time.*um_opt2(1))+(time.^2 .* um_opt2(2));

    r1(:,k)=sum(bestappr(:,k)-rmsrow);
    sigmasqr(k)=(1/length(rmsnew))*  r1(:,k)'*r1(:,k);
    r2(:,k)=sum(bestappr2(:,k)-rmsrow);
    sigmasqr2(k)=(1/length(rmsnew))*  r2(:,k)'*r2(:,k);

end

a=1;
for i=1:1
    axis=dt*linspace(1,length(rmsrow),length(rmsrow));
    figure('units','normalized','outerposition',[0 0 1 1])
    num=1:numel(rmsnew(1:length(rmsrow),a));
    hold on
    plot(axis(1:1:end),rmsnew(1:1:end,a),'.','MarkerSize',10)
    plot(axis,bestappr(1:length(rmsrow),a),'LineWidth',2)
    plot(axis,bestappr2(1:length(rmsrow),a),'LineWidth',2)
    % plot(rms)
    legend('Data','Linear Fit','Quadratic Fit')
    title('Sum of Root Mean Squares For Diffusing Particles')
    xlabel('Time(s)')
    ylabel('Sum of Distance Squared(\mu m^2)')
    h = gca;
    set(h,'FontSize',[32]);
    hold off
    a=a+1;
end


%% Plotting of raw data
a=1;
for i=1:7
    axis=dt*linspace(1,beancount(a)-1,beancount(a)-1);
    figure('units','normalized','outerposition',[0 0 1 1])

    hold on
    plot(axis(3:3:end),rms(3:3:beancount(a)-1,a),'.','MarkerSize',10)
    plot(axis,bestapprox(1:beancount(a)-1,a),'LineWidth',1.8)
    plot(axis,bestapprox2(1:beancount(a)-1,a),'LineWidth',1.8)
    % plot(rms)
    legend('Diffusion Data','Linear Fit','Quadratic Fit')
    title('Sum of Root Mean Squares For Diffusing Particles')
    xlabel('Time(s)')
    ylabel('Sum of Distance Squared(\mu m^2)')
    h = gca;
    set(h,'FontSize',[24]);
    hold off
    a=a+1;
end
%% 

a=1;
rs(:,a)=bestapprox2(:,a)-rms(:,a);
figure('units','normalized','outerposition',[0 0 1 1])

hold on
axis=dt*linspace(1,beancount(a)-1,beancount(a)-1);
plot(axis,rs(1:beancount(a)-1,a),'.','LineWidth',1.8)
% plot(rms)
legend('Data','Linear Fit','Exponential Fit')
title('X position fit')
xlabel('Time(s)')
ylabel('Sum of Distance Squared(um)')
h = gca;
set(h,'FontSize',[24]);
hold off

%% Functions


function J = linfun(um,obs,rho)

y=um(1)+(rho.*um(2));
J = (y-obs)'*((y-obs));

end 

function J = quadfun(um,obs,rho)

y=um(1)+(rho.*um(2))+(rho.^2 .*um(3)) ;

J = (y-obs)'*((y-obs));

end


function J = linfun2(um,obs,rho)

y=(rho.*um);
J = (y-obs)'*((y-obs));

end 

function J = quadfun2(um,obs,rho)

y=(rho.*um(1))+(rho.^2 .*um(2)) ;

J = (y-obs)'*((y-obs));

end





