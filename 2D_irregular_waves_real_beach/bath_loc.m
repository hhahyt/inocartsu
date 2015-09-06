clear all
close all
clc


%% IMPORTANT - decision variable

 % If loading data from already created data files,
%% then set "load_type" equal to 1, if you are going to create
%% a new bathymetry using this script, then set "load_type" to 2

%*****************************************************************
%%%%%%%%%%%%%%%   THIS IS THE SECTION WHERE THE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   2D BATHYMETRY IS CREATED      %%%%%%%%%%%%%%%%%%


%% If loading an already created data file, load here:
load nested_grid.grd
h=nested_grid;  % alter this line so that "h", the water

d=13; % Local depth of the wavemaker
p=13; %(Xmin – 510) between (Xmin – 120).
q=18; %Create an area to create a wavemaker.(p*dx=17*30=510)

[m,n]=size(h);

type=input('Enter 1 to use mean sea level(MSL)or\nEnter 2 to superimpose storm surge waves over MSL(100year mean = 2.6m)\n','s');
   
%% Normalize the values to mean sea level(MSL) 
%(MHW-MSL = 0.132) - MHW=mean high water
for j=1:m;
     for i=1:n;
       MSL(j,i)=h(j,i)+0.132;
     end
end
%% Superimpose storm surge waves over MSL (100year mean = 2.6m)
if type == '2';
    for j=1:m;
        for i=1:n;
         SS(j,i)=MSL(j,i)-2.6;
        end
    end
end
%********************************************************************
if type == '2';
  H=SS;
else
  H=MSL;
end
%% Interpolation boundaries
x0=390.00;
x1=510.00;

% 
%% Specify the starting point x_min and y_min

Xmin=0;  % Subtract 510m (17 grids) to create the wavemaker
Xmax=4620.00;  % Multiply dx*129=3840 is distance in meters of x in the grid.

Ymin=0;
Ymax=4530.00; % Multiply dx*72=2130 is distance in meters of y in the grid.

%% Specify grid lengths:
dx=30.00;
dy=30.00;

% Prepare matrices of new coordinates at which to interpolate

   
    [XI,YI] = meshgrid(Xmin:dx:Xmax,Ymin:dy:Ymax);

[mm,nn]=size(XI);

%% Extend the subset to create the an area for the wavemaker

for r=1:mm;
  for c=1:nn;
    if c <= p
        depth(r,c)=-d;
    elseif ((c > p) & (c < q))
        depth(r,c)=NaN;
    elseif c >= q
        depth(r,c)= H(r,c-17);
    end
  end
end

%% Interpolation of 120 meters(4 grids) between constant depth and real
%% bathymetry

disp('Computing an extended area for the wavemaker')

C13=depth(:,p);
C18=depth(:,q);

A=[C13 C18];

[wx,wy]=meshgrid((x0-dx):150:x1,Ymin:dy:Ymax);

[wxi,wyi]=meshgrid((x0-dx):dx:x1,Ymin:dy:Ymax);

wzi=interp2(wx,wy,A,wxi,wyi); %intepolation the area between 

cc=1;
for c=p:q;
    depth(:,c)=wzi(:,cc);
    cc=cc+1; 
end

%% Creating the grid

for j=1:mm;
   y(j)=Ymin+((j-1)*dy);
end

for i=1:nn;
    x(i)=Xmin+((i-1)*dx);
end

   
%% plot the data

figure(1)
[C,w]=contourf(x,y,depth,'LevelStep',2);
clabel(C,w)
hold on
v=[0 0];
% [C,w]=contour(x,y,depth,v,'k-','LineWidth',2);
% hold on
a1=167691;   b1=42570;    
a2=167661;   b2=42600;  
a3=167631;   b3=42630;  
a4=167721;   b4=42630;    
a5=167601;   b5=42750; 
a6=167451;   b6=43020; 
plot(a1,b1,'w*',a2,b2,'w*',a3,b3,'w*',a4,b4,'w*',a5,b5,'w*',a6,b6,'w*');
hold on
view(0,90)
grid off
colorbar
title('Playa Buye, Cabo Rojo')
xlabel('EASTINGS(m)')
ylabel('NORTHINGS (m)')

%*******************************************************************


figure(2)
surf(wxi,wyi,wzi)
hold on
% view(-75,50)
grid off
colorbar
%% Transform and save the data so that it can be loaded by the numerical program. The saves files are
%% called x_topo.dat, y_topo.dat, f_topo.dat, and size_topo.dat, and must
%% be moved to the same diretory as the numerical simulation.


x=x';
y=y';


save x_topo.dat x -ascii
save y_topo.dat y -ascii

sx=[size(x,1),size(y,1)]';

save size_topo.dat sx -ascii

m=length(x);
n=length(y);

file_dat=reshape(depth',[n*m,1]);

save f_topo.dat file_dat -ascii

% These saved to disk output files need to be moved to
% the directory from which the exectuable is run
