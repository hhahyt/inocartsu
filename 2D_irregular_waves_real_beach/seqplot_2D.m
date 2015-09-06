clear all

load time.dat
nt=length(time);
if nt==0
   nt=10000;
end

load s.dat

N_levels=s(2);
num_procs=s(3);

for rank=0:num_procs-1
   ind=['00' num2str(rank)];
   file_ind=ind(length(ind)-2:length(ind));
   
   eval(['load xloc' file_ind '.dat'])
   eval(['load yloc' file_ind '.dat'])
   
   filename=['dpth' file_ind '.dat'];
   fid(5*rank+1)=fopen(filename, 'r');
   
   filename=['zeta' file_ind '.dat'];
   fid(5*rank+2)=fopen(filename, 'r');         % open the file
   
   filename=['velo' file_ind '.dat'];
   fid(5*rank+3)=fopen(filename, 'r');         % open the file
   
   filename=['blvs' file_ind '.dat'];
   fid(5*rank+4)=fopen(filename, 'r');        % open the file  
   
   filename=['vort' file_ind '.dat'];
   fid(5*rank+5)=fopen(filename, 'r');        % open the file
end



% plot initial depth and determine glocal axis limits
clf
mx_x=0;
mn_x=1e6;
mx_y=0;
mn_y=1e6;

for rank=0:num_procs-1
   ind=['00' num2str(rank)];
   file_ind=ind(length(ind)-2:length(ind));
   
   eval(['x=xloc' file_ind ';'])
   eval(['y=yloc' file_ind ';'])
   
   mx_x=max(mx_x,max(x));
   mx_y=max(mx_y,max(y));
   
   mn_x=min(mn_x,min(x));
   mn_y=min(mn_y,min(y));
   
   nx=length(x);
   ny=length(y);
   
   ho=zeros(nx,ny);
   
   dum=fread(fid(5*rank+1),1,'int32');
   ho(:,:)=fread(fid(5*rank+1),[nx,ny],'single');
   dum=fread(fid(5*rank+1),1,'int32');
   
   surf(x,y,transpose(ho(:,:)))
   view(0,90)
   shading flat
   
   hold on
end
colorbar
pause

% plot various surfaces
for n=1:100
   clf
   [n,nt]
   for rank=0:num_procs-1
      ind=['00' num2str(rank)];
      file_ind=ind(length(ind)-2:length(ind));
      
      eval(['x=xloc' file_ind ';'])
      eval(['y=yloc' file_ind ';'])
      
      nx=length(x);
      ny=length(y);
      
      zeta=zeros(nx,ny,N_levels);
      tmp=zeros(nx,ny);
      bl_hor_wall=zeros(nx,ny);
      t_break=zeros(nx,ny);
      h=zeros(nx,ny);
      u=zeros(nx,ny,N_levels);
      v=zeros(nx,ny,N_levels);
      vort=zeros(nx,ny,N_levels);
      us=zeros(nx,ny,N_levels);
      vs=zeros(nx,ny,N_levels);
      
      for s=1:N_levels
%         dum=fread(fid(5*rank+3),1,'int32');
%         u(:,:,s)=fread(fid(5*rank+3),[nx,ny],'single');   % read x-velocity array, at z_alpha
%         v(:,:,s)=fread(fid(5*rank+3),[nx,ny],'single');   % read y-velocity array, at z_alpha
%         dum=fread(fid(5*rank+3),1,'int32');
         
         dum=fread(fid(5*rank+2),1,'int32');
         zeta(:,:,s)=fread(fid(5*rank+2),[nx,ny],'single');   % read free surface array
         dum=fread(fid(5*rank+2),1,'int32');         
         
 %        dum=fread(fid(5*rank+5),1,'int32');
 %        vort(:,:,s)=fread(fid(5*rank+5),[nx,ny],'single');   % read surface vorticity array
 %        us(:,:,s)=fread(fid(5*rank+5),[nx,ny],'single');   % read x-velocity array, at zeta
 %        vs(:,:,s)=fread(fid(5*rank+5),[nx,ny],'single');   % read y-velocity array, at zeta         
 %        dum=fread(fid(5*rank+5),1,'int32');       
      end 
      
      
      dum=fread(fid(5*rank+4),1,'int32');
      bl_hor_wall=fread(fid(5*rank+4),[nx,ny],'int32');  % location of wet=0 and dry=99 cells
      t_break=fread(fid(5*rank+4),[nx,ny],'single');  % breaking locations
      h=fread(fid(5*rank+4),[nx,ny],'single');  % water depth, only different from h0 for landslide
      dum=fread(fid(5*rank+4),1,'int32');
      
   %   subplot(2,1,1)
      v=[0:.5:5];
      contour3(x,y,h',v,'k')
      hold on
      z=transpose(zeta(:,:,1).*(1-bl_hor_wall/99)-(bl_hor_wall/99)*0);
      surf(x,y,z)
      
    %  subplot(2,1,2)
    %  plot(x,zeta(:,1,1),x,-h(:,1))
      
   end
   
     %    subplot(2,1,1)
   shading flat
   axis equal
   axis([mn_x mx_x mn_y mx_y -10 10])
   xlabel('x (m)')
   ylabel('y (m)')
   %title(['Time (sec) = ', num2str(time(n))])
   view(0,90)
%   caxis([-.5 .75])
   colorbar
   pause(.01)
end

fclose(fid);

