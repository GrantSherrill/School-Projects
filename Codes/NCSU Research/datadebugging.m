% % Main Block where data is read in, and processed to fit the RMS to a
% linear model and quadratic model. The coefficient of fit in the quadratic
% model is the velocity of the bead, which goes into stokes drag.
clear all;
close all;
filename='ResultsFromBead3mosaic.csv';

mu=300.006*10^(-3)*10^-12; %In N/m^2*s    *10^-12(m^2/um^2)
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

j=1;
k=1;


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
dim=size(beadx);
i=1;
for m=1:dim(2)
    if beadx(20,m)==0
        falsetrack(i)=m;
        i=i+1;
    end
end

for i=1:length(falsetrack)
    beadx(:,falsetrack(i)-i+1)=[];
    beady(:,falsetrack(i)-i+1)=[];
end
count=1;
newdim=size(beadx);
for n=1:newdim(2)
    for i=1:length(beadx)
        if  abs(beadx(i,n)-(beadx(1,n))) > 0.5 || abs(beady(i,n)-beady(1,n)) > 1
            goodtrack(count)=1;
            
        end
    
    end
    count=count+1;
end




% 
% off=30;
% 
% beadx(off,:)=[];
% for k=3:4
%     if beancount(k)<20
%         beadx(:,k)=0; 
%     end 
%     
%     
%     if beancount(k)<0
%        beancount(k)=1 ;
%     end
% end

beadxalg=beadx;
for k=1:newdim(2)
      tempbeadx=beadx(1:newdim(2)-1,k);
      tempbeady=beady(1:newdim(2)-1,k);
      tempdisp=sqrt(tempbeadx.^2+tempbeady.^2);
      vars(k)=var(tempbeadx);



  
    for j=1:newdim(2)-1
        newpos(j,k)=(sqrt((tempbeadx(j)-tempbeadx(1))^2 + (tempbeady(j)-tempbeady(1))^2))^2;
        rms(j,k)=sum(newpos(1:j,k))/j;
    end
    rmstemp=rms(1:newdim(2)-1,k);

    umguess=[1];
    time=dt*linspace(1,newdim(2)-1,newdim(2)-1);
    time=time';
    modelfun = @(um) linfun2(um , rmstemp,time);
    [um_opt] = fminsearch(modelfun,umguess);
    slope(k)=um_opt(1);

    
    bestapprox(1:length(rmstemp),k) = (time.*um_opt(1)); %Linear fit



    umguess2=[1,1];
    time=dt*linspace(1,newdim(2)-1,newdim(2)-1);
    time=time';
    modelfun2 = @(um2) quadfun2(um2 , rmstemp,time);
    [um_opt2] = fminsearch(modelfun2,umguess2); %quadratic fit (ax+b^2 x^2)
    slope2(k)=um_opt2(1);
    diffforce(k)=6*pi*mu*r*sqrt((um_opt2(2))); %Forces reported in model





    if um_opt2(2)<0   %Prevents imaginary forces in the case of plotting forces
        diffforce(k)=-imag(diffforce(k));
    end    
    bestapprox2(1:length(rmstemp),k) = (time.*um_opt2(1))+(time.^2 .* um_opt2(2));

%     Error Caluclations
    
    

    x=[rmstemp];
    rs=bestapprox(1:newdim(2)-1,k)-rmstemp;
    sigmasqr(k)=(1/(newdim(2)-2))*  rs(:)'*rs(:);
    V(k)=sigmasqr(k)*inv(x'*x);
    liner=sqrt(V(1,1));

end




%% 
a=1;
for i=1:20
    axis=dt*linspace(1,newdim(2)-1,newdim(2)-1);
    figure('units','normalized','outerposition',[0 0 1 1])

    hold on
    plot(axis(1:1:end),rms(1:1:newdim(2)-1,a),'.','MarkerSize',10)
    plot(axis,bestapprox(1:newdim(2)-1,a),'LineWidth',1.8)
    plot(axis,bestapprox2(1:newdim(2)-1,a),'LineWidth',1.8)
    % plot(rms)
    legend('Diffusion Data','Linear Fit','Quadratic Fit')
    title('Sum of Root Mean Squares ')
    xlabel('Time(s)')
    ylabel('Sum of Distance Squared(\mu m^2)')
    h = gca;
    set(h,'FontSize',[24]);
    hold off
    a=a+1;
end



%% 


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





