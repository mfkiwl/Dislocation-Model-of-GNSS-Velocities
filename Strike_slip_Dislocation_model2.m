% The matlab code help to model the observed velocities from GPS stations across strike-slip faults 
% by analyzing misfit between the fault Slip rate, Locking Depth, Fault Dip, and Vertical offset (Based on methods by Okada 1992)

% Misfit analysis between the observed and modelled GPS velocities for strike-slip faults 
% (Variable: Slip rate, Locking Depth, Fault Dip & Vertical Offset uY)

% Last modified on: 17 Mar, 2026 by D. Panda

clear all
close all
clc
addpath(genpath(pwd))

%%

% Misfit analysis between observed and modelled GPS velocities for dip-slip faults
% (Variable: Slip rate, Locking Depth, Dip, Vertical Offset uY)

[file,location] = uigetfile('*');

Fault_slip=readtable(fullfile(location,file));

if width(Fault_slip)<3
    error('File should contain at least 3 columns (Distance, Displacement, Error)')
else
    disp("Input File contains 3 columns > Distance, Displacement, Error")

    %% Define the variables        
Fault_slip=sortrows(Fault_slip);
dist=table2array(Fault_slip(:,1))*(-1);
% dist=table2array(Fault_slip(:,1));
slip=table2array(Fault_slip(:,2));
error=table2array(Fault_slip(:,5));
plot(dist,slip,'-bdiamond','linewidth',1,'DisplayName','Observed')
hold on
errorbar(dist,slip,error,'.b','DisplayName','Observed')
title('Observed Displacement (mm/year)')
% text(gps_nrsc.Long,gps_nrsc.Lat,gps_nrsc.Station,'Vert','bottom', 'Horiz','left','FontSize',7)
end
%%

dist_new = dist;
slip_new = slip;
error_new = error;

menuop = menu('EXCLUSION OF ONE OR MORE POINTS FROM CALCULATIONS','YES','NO');

if menuop ==1
    
    figure(1)
    plot(dist,slip,'-bdiamond','linewidth',1,'DisplayName','Observed')
    hold on
    errorbar(dist,slip,error,'.b','DisplayName','Observed')

    [rx,ry]=getpts(figure(1));
    dist_range=2;
    slip_range=2;
    
    for k = 1:length(rx)
        id = [];
        id2 = [];
        id = find(dist<=fix(rx(k))+dist_range & dist>=fix(rx(k))-dist_range);
        id2 = find(slip<=fix(ry(k))+slip_range & slip>=fix(ry(k))-slip_range);
        if k > 1
            id = id - (k-1);
        end
        
        if length(id)<length(id2)
            
            plot(dist,slip,'-bdiamond','linewidth',1,'DisplayName','Observed')
            hold on
            plot(dist(id),slip(id),'-rx','linewidth',1,'DisplayName','Observed')

        dist_new(id,:)=[];

        slip_new(id,:)=[];

        error_new(id,:)=[];

    else

        dist_new(id2,:)=[];

        slip_new(id2,:)=[];

        error_new(id2,:)=[];

    figure (2)
    plot(dist,slip,'-bdiamond','linewidth',1,'DisplayName','Observed')
    hold on
    plot(dist(id2),slip(id2),'-rx','linewidth',1,'DisplayName','Observed')

    end
    end

else

    dist_new=dist;
    slip_new=slip;
    error_new=error;

end

hold off
figure (3)
plot(dist,slip,'-bdiamond','linewidth',1,'DisplayName','Actual Stations')
hold on
plot(dist_new,slip_new,'-rdiamond','linewidth',1,'DisplayName','After Removal')
legend

%%  Estimate the RMSE misfit

RMSE1 = [];
type='S';
Mu = 30e9;
Poisson = 0.25;
nobs = length(dist_new);

fprintf('Estimating RMSE by varying fault Slip rate, Locking Depth, Dip, Vertical Offset...\n')
tic
for i=0:500:40000     % Locking depth (in meters) 
    for j=0:1:25       % Fault slip rate (in mm)
        for k=1:45     % Fault Dip (in degrees)
            for m=-5:5 % Fault slip vertical offset (in mm)
            
fault=[-1e9,0;1e9,0];
dip1=k;
dip=(180-dip1)*pi/180;
depth=[0e3,i];
B=j;
y=dist_new*1e3;
x=0*y;
z=0*y;

[uX,uY,uZ]=Okada1992(x,y,z,fault,dip,depth,B,type,Mu,Poisson);
% uY(y>0)=uY(y>0)+(B*cosd(dip1));
uX(y<0) = uX(y<0)+(B);

% uY=uY+m;

% RMSE estimation between observed and modelled GPS velocities
%    RMSE1 = [RMSE1;i/1000,(i/1000)/sind(dip),j,sqrt(mean((slip(:)-uY(:)).^2))];
    RMSE1 = [RMSE1;i/1000,j,k,m,(1/nobs)*sqrt(sum((slip_new(:)-uX(:)).^2./error_new.^2))];

            end

            uX=[];

        end
    end
end
toc
fprintf('Estimating RMSE by varying fault Slip rate, Locking Depth, Dip, Vertical Offset...Done\n')

%%

x2=RMSE1(:,1);
y2=RMSE1(:,2);
z2=RMSE1(:,3);
p1=RMSE1(:,4);
p2=RMSE1(:,5);

minError = find(p2==min(p2));
depth2=x2(minError);
slip2=y2(minError);
dip2=z2(minError);
dip3=(180-dip2)*pi/180;
Vert_offset=p1(minError);

figure(1),clf
[uX,uY,uZ]=Okada1992(x,y,z,fault,dip3,[0,depth2*1e3],slip2,type,Mu,Poisson);
% uY(y>0)=uY(y>0)+(slip2*cosd(dip2));
uX(y<0) = uX(y<0)+(slip2);
plot(y/1e3,uX+Vert_offset,'-bx','linewidth',1,'DisplayName','Modelled')
hold on
errorbar(dist,slip,error,'DisplayName','Observed')
xlim([-300,300])
legend ('location','northwest')


figure(2),clf
nobs=200;
y=linspace(-500e3,500e3,nobs*2);
x=0*y;
z=0*y;
[uX,uY,uZ]=Okada1992(x,y,z,fault,dip3,[0,depth2*1e3],slip2,type,Mu,Poisson);
% uY(y>0)=uY(y>0)+(slip2*cosd(dip2));
uX(y<0) = uX(y<0)+(slip2);

legend(plot(y/1e3,uX+Vert_offset,'g','linewidth',1,'DisplayName','Best Fit'))
hold on
errorbar(dist,slip,error,'DisplayName','Observed')
xlim([-300,300])
legend ('location','northwest')
caption = sprintf('Dip = %g, Locking Depth = %g km, Slip = %g mm/yr, uY = %g mm/yr', dip2, depth2, slip2, Vert_offset);
title(caption, 'FontSize', 10);
