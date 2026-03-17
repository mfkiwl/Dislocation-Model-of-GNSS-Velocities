% The matlab code help to model the observed velocities from GPS stations across strike-slip faults 
% by analyzing misfit between the fault Slip rate & Locking Depth (Based on methods by Okada 1992)

% Misfit analysis between the observed and modelled GPS velocities for strike-slip faults
% (Variable: Slip rate & Locking Depth; Constant: Dip & Vertical Offset uY)

% Input: File must contain three columns (atleast) 
% Distance from Fault Trace, Displacement, Uncertainity

% Last modified on: 17 Mar, 2026 by D. Panda


clear all
close all
clc
addpath(genpath(pwd))

%%
% Load GNSS data
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
error=table2array(Fault_slip(:,3));
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

% Misfit between observed and modelled GNSS velocities for dip-slip faults
% (Variable: Slip rate & Locking Depth ; Constant: Dip & Vertical Offset uY)

plot(dist_new,slip_new,'-rdiamond','linewidth',1,'DisplayName','After Removal')
hold on
errorbar(dist_new,slip_new,error_new,'.r','DisplayName','Observed')
title('Final Observed Displacement (mm/year)')

% Modify hard coded parameters
RMSE1 = [];
type = 'S'; % Type of fault (Dip slip: 'D', Strike slip: 'S')
dip1 = 85;  % Dip of the fault
Mu = 30e9;  % Shear Modulus
Poisson = 0.25;
nobs = length(dist_new);    % Number of GPS stations
LDival = 100;     %    Locking Depth interval, in meters
slipival = 0.2;    %   Slip rate interval, in mm

fprintf('Estimating RMSE by varying fault slip rate and locking depth...\n')
tic
for i = 0:LDival:60000       % Locking depth (in meters)
    
    for j = 0:slipival:40     % Fault slip rate (in mm)
        
        fault = [-1e9,0;1e9,0];
        dip = (180-dip1)*pi/180;
        depth = [0e3,i];
        B = j;
        y = dist_new*1e3;
        x = 0*y;
        z = 0*y;

        [uX,uY,uZ] = Okada1992(x,y,z,fault,dip,depth,B,type,Mu,Poisson);
%         uX(y>0) = uY(y>0)-(B*cosd(dip1));
        uX(y<0) = uX(y<0)+(B);

        % RMSE estimation between observed and modelled GPS velocities
        %    RMSE1 = [RMSE1;i/1000,(i/1000)/sind(dip),j,sqrt(mean((slip(:)-uY(:)).^2))];
        RMSE1 = [RMSE1;i/1000,j,(1/nobs)*sqrt(sum((slip_new(:)-uX(:)).^2./error_new.^2))];
    end

    uX=[];

end
toc
fprintf('Estimating RMSE by varying fault slip rate and locking depth...Done\n')

%%

x2=RMSE1(:,1);
y2=RMSE1(:,2);
z2=RMSE1(:,3);
xival=0.2;      % x interval for the mesh
yival=0.2;      % y interval for the mesh
[X,Y]=meshgrid(min(x2):xival:max(x2),min(y2):yival:max(y2));
Z=griddata(x2,y2,z2,X,Y);

figure(3),clf
contourf(X, Y, Z)
minError = find(Z == min(Z(:)));
minX=X(minError);
minY=Y(minError);
hold on
plot(minX,minY,'o','color','red')

figure(4),clf
[uX,uY,uZ]=Okada1992(x,y,z,fault,dip,[0,minX*1e3],minY,type,Mu,Poisson);
% uY(y>0)=uY(y>0)+(minY*cosd(dip1));
% uY=uY+1;
uX(y<0) = uX(y<0)+(minY);
plot(y/1e3,uX,'-bx','linewidth',1,'DisplayName','Modelled')
hold on
errorbar(dist_new,slip_new,error_new,'DisplayName','Observed')
xlim([-500,500])
legend ('location','northwest')

figure(5),clf
nobs=200;
y=linspace(-500e3,500e3,nobs*2);
x=0*y;
z=0*y;
[uX,uY,uZ]=Okada1992(x,y,z,fault,dip,[0,minX*1e3],minY,type,Mu,Poisson);
% uY(y>0)=uY(y>0)+(minY*cosd(dip1));
uX(y<0) = uX(y<0)+(minY);

legend(plot(y/1e3,uX,'g','linewidth',1,'DisplayName','Best Fit'))
hold on
errorbar(dist_new,slip_new,error_new,'DisplayName','Observed')
xlim([-300,300])
legend ('location','northwest')

title(['Dip =  ' num2str(dip1), ' ; Locking Depth (km) = ' num2str(minX), ' ; Slip (mm) = ' num2str(minY)])
% title(['Dip =  ' num2str(dip1), ' ; Locking Width (km) = ' num2str(minX/tand(dip1)), ' ; Slip (mm) = ' num2str(minY)])

%%
% Write RMSE and Okada Locking curve files
fprintf('Saving file: Misfit between fault slip rate and locking depth...Done\n')
% save('SS_RMSE.txt','RMSE1','-ascii');
header1={'Locking_Depth_km','Slip_mm','rmse'}; 
file1=[header1;num2cell(RMSE1)];
writecell(file1,'SS_RMSE.txt');

fprintf('Saving file: Slip rate deficit curve...Done\n')
% save('SS_Okada_curve.txt','[(x/1000)',uY']','-ascii');
header2={'Distance_km','Slip_mm'}; 
file2=[header2;num2cell([(x/1000)',uY'])];
writecell(file2,'SS_Okada_curve.txt');


