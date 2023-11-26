% % Main Block where data is read in, and processed to fit the RMS to a
% linear model and quadratic model. The coefficient of fit in the quadratic
% model is the velocity of the bead, which goes into stokes drag.
clear all;
close all;
filename='Results from MGSTriprobecontrol.csv';

mu=0.0009*10^-12; %In N/um^2*s    *10^-12(m^2/um^2)
r=1; %in um

fr=19.38; %frames per second

data=importdata(filename);
time=data.data(:,1);
bead=data.data(:,2);
xpos=data.data(:,4);
ypos=data.data(:,5);

xpos=xpos/8.259; %From pixels to microns
ypos=ypos/8.259;

time=time/fr;
dt=1/fr;

T=299;






cutoff=15/50;



for k=1:length(bead)
    j=1;
    for i=1:length(bead)
        if bead(i)==k
            beadx(j,k)=xpos(i);
            beady(j,k)=ypos(i);
            beancount(k)=j;
            j=j+1;
            
        end
    end
end


%%%% Test for removing small columns
j=1;
k=1;
dim=size(beadx);
i=1;
for m=1:dim(2)
    if beadx(50,m)==0
        falsetrack(i)=m;
        i=i+1;
    end
end

for i=1:length(falsetrack)
    beadx(:,falsetrack(i)-i+1)=[];
    beady(:,falsetrack(i)-i+1)=[];
    beancount(falsetrack(i)-i+1)=[];
end



for k=1:length(beancount)

      tempbeadx=beadx(1:beancount(k)-1,k);
      tempbeady=beady(1:beancount(k)-1,k);
      vars(k)=var(tempbeadx);



    for m=1:beancount(k)-1
        for j=1:beancount(k)-1
            newpos(j,k)=(sqrt((tempbeadx(j)-tempbeadx(1))^2 + (tempbeady(j)-tempbeady(1))^2))^2;
            temppos(j,m)=(sqrt((tempbeadx(j)-tempbeadx(m))^2 + (tempbeady(j)-tempbeady(m))^2))^2;
            rms(j,k)=sum(newpos(1:j,k))/j;
            
    
    
        end
    end
    for n=2:beancount(k)-1
        vec=diag(temppos,n-1);
        vec=vec(1:(n-1):end);
        ptrms(n-1,k)=sum(vec)/length(vec);

        %ptrms(n,k)=sum(diag(temppos,n-1))/(beancount(k)-n+1);
%         testvec=diag(temppos,1);
%         vecslice=testvec(1:2:end);
        ptest(n-1,k)=std(vec)/sqrt(length(vec));
        ptest2(n-1,k)=sqrt(sum((vec-mean(vec)).^2)/length(vec));
    end
     ptrms(1:beancount(k)-2,k)=ptrms(1:beancount(k)-2,k)-ptrms(1,k);
    rmstemp=rms(1:beancount(k)-2,k);
    rmstemp2=ptrms(1:beancount(k)-2,k);



    umguess=[1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    short(k)=round(cutoff*(beancount(k)-1));
    sttime=dt*linspace(1,short(k),short(k));
    sttime=sttime';
    rmsshort=rmstemp2(1:short(k));
    modelfun = @(um) linfun2(um , rmsshort,sttime);
    [um_opt] = fminsearch(modelfun,umguess);
    slope(k)=um_opt(1);
    D(k)=um_opt(1) / 4; %in um^2/s 
    kb=1.38*10^(-11); %in um^2 kg s^-2 K-1
  
    
  

    Dest=kb*T/(6*pi*r*mu)*10^-6;
    
    visc(k)=kb*T / (6*pi*r*D(k))*10^6; % in kg  /um s, multiply by 10^6 um/m to get into Ns/m^2   
    radius(k)=kb*T / (6*pi*mu*D(k))*10^-6;

    

    bestapprox(1:length(rmsshort),k) =(sttime.*um_opt(1));

    umguess2=[1,1];
    modelfun2 = @(um2) quadfun2(um2 , rmsshort,sttime);
    [um_opt2] = fminsearch(modelfun2,umguess2);
    slope2(k)=um_opt2(1);
    quadslope(k)=um_opt2(2);
    diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));
    if um_opt2(2)<0
        diffforce(k)=-imag(diffforce(k));
    end    
    diffforce(k)=round(diffforce(k),3,"significant")*10^12;
    bestapprox2(1:length(rmsshort),k) = (sttime.*um_opt2(1))+(sttime.^2 .* um_opt2(2));

%Error
%     Error Caluclations
%     
    
    x2=[sttime];
    rs2(1:length(rmsshort),k)=bestapprox2(1:short(k),k)-rmsshort;
    sigmasqr2(k)=(1/(short(k)-2))*  rs2(1:length(rmsshort),k)'*rs2(1:length(rmsshort),k);
    terror2=sigmasqr2(k).*inv(x2'*x2);
    liner(k)=sqrt(terror2(1,1));
 


    x=[sttime,sttime.^2];
    rs(1:length(rmsshort),k)=bestapprox2(1:short(k),k)-rmsshort;
    sigmasqr(k)=(1/(short(k)-2))*  rs(1:length(rmsshort),k)'*rs(1:length(rmsshort),k);
    terror=sigmasqr(k).*inv(x'*x);
    par1er(k)=sqrt(terror(1,1));
    par2er(k)=sqrt(terror(2,2));
    differ(k)=6*pi*mu*r*sqrt((par2er(k)))*10^12;


clear('tempbeadx','tempbeady','temppos')
end

%% plotting
a=1;
for i=1:4
    axis=dt*linspace(1,beancount(a)-2,beancount(a)-2);
    staxis=dt*linspace(1,short(a),short(a));
    figure('units','normalized','outerposition',[0 0 1 1])

    hold on
    plot(axis(1:end),ptrms(1:beancount(a)-2,a),'.','MarkerSize',10)
    plot(staxis,bestapprox(1:short(a),a),'LineWidth',1.8)
    plot(staxis,bestapprox2(1:short(a),a),'LineWidth',1.8)
    errorbar(axis,ptrms(1:1:beancount(a)-2,a),ptest(1:1:beancount(a)-2,a),'k','LineStyle','none')
    % plot(rms)
    legend('Diffusion Data','Linear Fit','Quadratic Fit')
    title('Sum of Root Mean Squares For 8mT')
    xlabel('\delta t(s)')
    ylabel('Sum of Distance Squared(\mu m^2)')
    h = gca;
    set(h,'FontSize',[24]);
    hold off

%     figure('units','normalized','outerposition',[0 0 1 1])
%     hold on
%     plot(axis(1:end),0*ptrms(1:beancount(a)-1,a),'LineWidth',1.6)
%     plot(staxis,rs(1:short(a),a),'.','MarkerSize',10)
%     hold off
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

%% Histogram of Forces
for k=1:length(beancount)
    if beancount(k)<20
        diffforce(k)=0;
    end
end
figure()
histogram(diffforce*10^15,50)
legend('Data')
title('X position fit')
xlabel('Force(fN)')
ylabel('Sum of Distance Squared(um)')
h = gca;
set(h,'FontSize',[24]);
hold off


%% Experimental Stuff

k=64;
s=100;
for j=1:s
        newpos(j,k)=((beadx(j,k)-beadx(1,k)))^2;
        rmsch(j,k)=sum(newpos(1:j,k))/j;
end
umguess2=[1,1];
time=dt*linspace(1,s,s);
time=time';
modelfun2 = @(um2) quadfun2(um2 , rmsch(:,k),time);
[um_opt2] = fminsearch(modelfun2,umguess2);
slope2(k)=um_opt2(1);
diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));


%% Bead Starting Position and Force origin videos

BOI=beady(:,1);
tb=dt*linspace(1,length(BOI),length(BOI));

starts=20;
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

a=3;
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


% 
% figure()
% hold on
% plot(tb,BOI,'.','LineWidth',1.8)
% % plot(rms)
% legend('Data')
% title('X position')
% xlabel('Time(s)')
% ylabel('Sum of Distance Squared(um)')
% h = gca;
% set(h,'FontSize',[24]);
% hold off


%% Y Component


for k=1:length(beancount)
      tempbeadx=beadx(1:beancount(k)-1,k);
      tempbeady=beady(1:beancount(k)-1,k);
      tempdisp=sqrt(tempbeadx.^2+tempbeady.^2);





  
    for j=1:beancount(k)-1
        newpos(j,k)=(tempbeady(j)-tempbeady(1))^2;
        rms(j,k)=sum(newpos(1:j,k))/j;
    end
    rmstemp=rms(1:beancount(k)-1,k);

    umguess=[1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    modelfun = @(um) linfun2(um , rmstemp,time);
    [um_opt] = fminsearch(modelfun,umguess);
    slope(k)=um_opt(1);

    
    bestapprox(1:length(rmstemp),k) = (time.*um_opt(1));

    umguess2=[1,1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    modelfun2 = @(um2) quadfun2(um2 , rmstemp,time);
    [um_opt2] = fminsearch(modelfun2,umguess2);
    slope2(k)=um_opt2(1);
    diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));
    if um_opt2(2)<0
        diffforce(k)=-imag(diffforce(k));
    end    
    bestapprox2(1:length(rmstemp),k) = (time.*um_opt2(1))+(time.^2 .* um_opt2(2));

%     Error Caluclations
    
    

%     x=[rmstemp,rmstemp.^2];
%     rs2=bestapprox2(1:beancount(k)-1,k)-rmstemp;
%     sigmasqr2=(1/(beancount(k)-2))*  rs2(:)'*rs2(:);
%     V=sigmasqr2*inv(x'*x);
%     liner=sqrt(V(1,1));
%     exper=sqrt(V(2,2));

end
%% Total Displacement

for k=1:length(beancount)
      tempbeadx=beadx(1:beancount(k)-1,k);
      tempbeady=beady(1:beancount(k)-1,k);
     





  
    for j=1:beancount(k)-1
       newpos(j,k)=(sqrt((tempbeadx(j)-tempbeadx(1))^2 + (tempbeady(j)-tempbeady(1))^2))^2;
        rms(j,k)=sum(newpos(1:j,k))/j;
    end
    rmstemp=rms(1:beancount(k)-1,k);

    umguess=[1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    modelfun = @(um) linfun2(um , rmstemp,time);
    [um_opt] = fminsearch(modelfun,umguess);
    slope(k)=um_opt(1);

    
    bestapprox(1:length(rmstemp),k) = (time.*um_opt(1));

    umguess2=[1,1];
    time=dt*linspace(1,beancount(k)-1,beancount(k)-1);
    time=time';
    modelfun2 = @(um2) quadfun2(um2 , rmstemp,time);
    [um_opt2] = fminsearch(modelfun2,umguess2);
    slope2(k)=um_opt2(1);
    diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2)));
    if um_opt2(2)<0
        diffforce(k)=-imag(diffforce(k));
    end    
    bestapprox2(1:length(rmstemp),k) = (time.*um_opt2(1))+(time.^2 .* um_opt2(2));
end



%% 
% beadxn=beadx;
% beadxn(342,:)=[];
a=1;
for i=1:2
    trackdiff(:,i)=beadxalg(:,i)-beadxn(1:length(beadxalg),(3-i));
    avgtrack(i)=mean(trackdiff(i));
end
for i=1:2
    axis=linspace(1,length(beadxalg(:,i)),length(beadxalg(:,i)));
    figure('units','normalized','outerposition',[0 0 1 1])

    hold on
    plot(axis,beadxalg(:,a),'LineWidth',1.8)
    plot(axis,beadxn(1:length(beadxalg),(3-a)),'LineWidth',1.8)
    % plot(rms)
    legend('Algorithm','Manual')
    title('X position tracking difference')
    xlabel('Time(s)')
    ylabel('Distance(um)')
    h = gca;
    set(h,'FontSize',[24]);
    hold off
    a=a+1;
end

%% 
a=1;
for i=1:2
    trackdiff(:,i)=beadxalg(:,i)-beadxn(1:length(beadxalg),(3-i));
    avgtrack(i)=mean(trackdiff(i));
end
for i=1:2
    axis=linspace(1,length(beadxalg(:,i)),length(beadxalg(:,i)));
    figure('units','normalized','outerposition',[0 0 1 1])

    hold on
    plot(axis,trackdiff(:,a),'.','MarkerSize',24)
    plot(axis,0*trackdiff(:,a),'LineWidth',1.8)
    plot(axis,0*trackdiff(:,a)+avgtrack(i),'LineWidth',1.8)
    % plot(rms)
    legend('Data','Zero','Mean')
    title('X position tracking difference')
    xlabel('Time(s)')
    ylabel('Distance(um)')
    h = gca;
    set(h,'FontSize',[24]);
    hold off
    a=a+1;
end


%% 


function J = linfun(um,obs,rho)

y=um(2)+(rho.*um(1));
J = (y-obs)'*((y-obs));

end 

function J = quadfun(um,obs,rho)

y=um(3)+(rho.*um(1))+(rho.^2 .*um(2)) ;

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





