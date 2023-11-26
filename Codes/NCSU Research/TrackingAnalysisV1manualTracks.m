%% Manual Tracking inputs
clear all
filename='Results from MGSTriprobe3vid2.csv';

mu=0.0009*10^-12; %In N/m^2*s    *10^-12(m^2/um^2)
r=1; %in um

fr=19.38; %frames per second

data=importdata(filename);
time=data.data(:,2);
bead=data.data(:,1);
xpos=data.data(:,3);
ypos=data.data(:,4);

xpos=xpos/8.26; %From pixels to microns
ypos=ypos/8.26;

time=time/fr;
dt=1/fr;


cutoff=40/100;

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

% 
% %%%% Test for removing small columns
% j=1;
% k=1;
% dim=size(beadx);
% i=1;
% for m=1:dim(2)
%     if beadx(20,m)==0
%         falsetrack(i)=m;
%         i=i+1;
%     end
% end
% 
% for i=1:length(falsetrack)
%     beadx(:,falsetrack(i)-i+1)=[];
%     beady(:,falsetrack(i)-i+1)=[];
%     beancount(falsetrack(i)-i+1)=[];
% end



for k=1:length(beancount)
      tempbeadx=beadx(1:beancount(k)-1,k);
      tempbeady=beady(1:beancount(k)-1,k);
      tempdisp=sqrt(tempbeadx.^2+tempbeady.^2);





  

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
        ptrms(n,k)=sum(vec)/length(vec);
        %ptrms(n,k)=sum(diag(temppos,n-1))/(beancount(k)-n+1);
%         testvec=diag(temppos,1);
%         vecslice=testvec(1:2:end);
        ptest(n,k)=std(vec)/sqrt(length(vec));
        ptest2(n,k)=sqrt(sum((vec-mean(vec)).^2)/length(vec));
   end
   ptrms(1:beancount(k)-2,k)=ptrms(1:beancount(k)-2,k)-ptrms(1,k);
    rmstemp=rms(1:beancount(k)-1,k);
    rmstemp2=ptrms(1:beancount(k)-1,k);


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
    T=300;
    
    visc(k)=kb*T / (6*pi*r*D(k))*10^6; % in kg  /um s, multiply by 10^6 um/m to get into Ns/m^2   
    radius(k)=kb*T / (6*pi*mu*D(k))*10^-6;

    bestapprox(1:length(rmsshort),k) = (sttime.*um_opt(1));

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


%%Error
   %%% Error Caluclations
    
    
    x=[sttime,sttime.^2];
    rs(1:length(rmsshort),k)=bestapprox2(1:short(k),k)-rmsshort;
    sigmasqr(k)=(1/(short(k)-2))*  rs(1:length(rmsshort),k)'*rs(1:length(rmsshort),k);
    terror=sigmasqr(k).*inv(x'*x);
    par1er(k)=sqrt(terror(1,1));
    par2er(k)=sqrt(terror(2,2));
    differ(k)=6*pi*mu*r*sqrt((par2er(k)))*10^12;

   

clear('tempbeadx','tempbeady','temppos')
end


%% Plotting of raw data
a=6;
for i=1:2
    axis=dt*linspace(1,beancount(a)-1,beancount(a)-1);
    staxis=dt*linspace(1,short(a),short(a));
    figure('units','normalized','outerposition',[0 0 1 1])

    hold on
    plot(axis(1:end),ptrms(1:beancount(a)-1,a),'.','MarkerSize',10)
    plot(staxis,bestapprox(1:short(a),a),'LineWidth',1.8)
    plot(staxis,bestapprox2(1:short(a),a),'LineWidth',1.8)
    errorbar(axis,ptrms(1:1:beancount(a)-1,a),ptest2(1:1:beancount(a)-1,a),'k','LineStyle','none')
    % plot(rms)
    legend('Diffusion Data','Linear Fit','Quadratic Fit')
    title('Sum of Root Mean Squares For Diffusing Particles')
    xlabel('\delta t(s)')
    ylabel('Sum of Distance Squared(\mu m^2)')
    h = gca;
    set(h,'FontSize',[24]);
    hold off
    a=a+1;
end
%% Histograms

b=sort(diffforce);

% nbins = 4;
% [N,edges]=histcounts(diffforce);
% h = histogram(diffforce,nbins);
% dpb=length(diffforce)/nbins;
for i=1:nbins
    er(i)=sum(differ(1+dpb*(i-1):dpb+dpb*(i-1)));
end
hold on
errorbar(diffforce,differ,'ko')

hold off
    

%% %NonBiased diffusion velocity

for   k=1:length(beancount) 
    for m=1:beancount(k)-2
        velx1(m,k)= fr* (beadx(m+1,k)-beadx(m,k))/2; %In microns per second
        vely1(m,k)= fr* (beady(m+1,k)-beady(m,k))/2; %In microns per second     
        veldis(m,k)=sqrt(velx1(m,k).^2+vely1(m,k).^2);
        
    end 
    
    velavg(k)=mean(veldis(:,k));
    maxvx=max(abs(velx1));
    minvx=min(abs(velx1));
    avgvx=mean((velx1));
    
    
    maxvy=max(vely1);
    minvy=min(vely1);
    avgvy=mean(vely1);
    
%     velx2=(beadx(length(beadx(:,k)))-beadx(1,k))/beancount(k); %In microns per second
%     vely2=(beady(length(beady(:,k)))-beady(1,k))/beancount(k); %In microns per second

    
    forcemax(k)=6*pi*r*mu*maxvy(k);
    forceavg(k)=6*pi*r*mu*avgvy(k);
%     forceavg2(k)=6*pi*r*mu*vely2;
    forcedispl(k)=6*pi*r*mu*velavg(k);

end
%% y displacement

clear all
filename='ResultsFromBead4.csv';

mu=330.06*10^(-3)*10^-12; %In N/m^2*s    *10^-12(m^2/um^2)
r=1; %in um

fr=19.38; %frames per second

data=importdata(filename);
time=data.data(:,2);
bead=data.data(:,1);
xpos=data.data(:,3);
ypos=data.data(:,4);

xpos=xpos/8.259; %From pixels to microns
ypos=ypos/8.259;

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
        newpos(j,k)=(sqrt( (tempbeady(j)-tempbeady(1))^2))^2;
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

        %NonBiased diffusion velocity

    
    for m=1:beancount(k)-2
    velx1(m)= (tempbeadx(m+1)-tempbeadx(m))*fr; %In microns per second
    vely1(m)= (tempbeady(m+1)-tempbeady(m))*fr; %In microns per second
    end 
    veldisp= sqrt(velx1.^2+vely1.^2 );
    
    maxvx=max(abs(velx1));
    minvx=min(abs(velx1));
    avgvx=mean((velx1));
    
    
    maxvy=max(vely1);
    minvy=min(vely1);
    avgvy=mean(vely1);
    

    avgv=mean(veldisp);

    velx2(k)=(tempbeadx(length(tempbeadx))-tempbeadx(1))/beancount(k); %In microns per second
    vely2(k)=(tempbeady(length(tempbeady))-tempbeady(1))/beancount(k); %In microns per second
    veldis(k)=sqrt(velx2(k)^2+vely2(k)^2);
    
    forcemax(k)=6*pi*r*mu*maxvx;
    forceavg(k)=6*pi*r*mu*avgv;
    forceavg2(k)=6*pi*r*mu*velx2(k);
    forcedispl(k)=6*pi*r*mu*veldis(k);


end
%% x data

clear all
filename='Results from 90glyc BigProbe8A319hzFinalTest_BrightfieldConfocal_2022-04-01_WidefieldWhiteLightBeadTimelapse_14.csv';

mu=234.6*10^(-3)*10^-12; %In N/m^2*s    *10^-12(m^2/um^2)
r=1; %in um

fr=19.38; %frames per second

data=importdata(filename);
time=data.data(:,2);
bead=data.data(:,1);
xpos=data.data(:,3);
ypos=data.data(:,4);

xpos=xpos/8.259; %From pixels to microns
ypos=ypos/8.259;

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

end
%% Bead Starting Position and Force origin videos

tempbeadx2=beadx(1:beancount(2)-1,2);
tempbeady2=beady(1:beancount(2)-1,2);
tb=dt*linspace(1,length(tempbeadx2),length(tempbeadx2));

starts=6;
off=1;

for m=1:length(tempbeadx2)
    if tempbeadx2(m)~=0
    ep=m;
    end
end
for k=1:starts
    for j=1:ep
        posm(j,k)=(sqrt((tempbeadx2(j)-tempbeadx2(k))^2 + (tempbeady2(j)-tempbeady2(k))^2))^2;
        rmsnew(j,k)=sum(posm(1:j,k))/j;
    end
end
% rmsnew(1:off,:)=[];

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
for i=1:6
    axis=dt*linspace(1,length(rmsrow),length(rmsrow));
    figure('units','normalized','outerposition',[0 0 1 1])
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



%% total displacement testimg

clear all
filename='ResultsFromBead1.csv';

mu=330.06*10^(-3)*10^-12; %In N/m^2*s    *10^-12(m^2/um^2)
r=1; %in um

fr=19.38; %frames per second

data=importdata(filename);
time=data.data(:,2);
bead=data.data(:,1);
xpos=data.data(:,3);
ypos=data.data(:,4);

xpos=xpos/8.259; %From pixels to microns
ypos=ypos/8.259;

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
        newpos(j,k)=(sqrt((tempbeadx(j)-tempbeadx(1))^2 + (tempbeady(j)-tempbeady(1))^2))^2;
        rms(j,k)=sum(newpos(1:j,k))/j;
    end
    for c=1:beancount(k)-1
        for j=1:beancount(k)-c
            newpos2(j,c)=(sqrt((tempbeadx(j+c-1)-tempbeadx(c))^2 + (tempbeady(j+c-1)-tempbeady(c))^2))^2;
            rms2(j,c)=sum(newpos2(1:j,c))/(j-c);
        end
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





