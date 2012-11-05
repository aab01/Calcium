function VisualiseCellSpikesV1(Cells)
% function VisualiseCellSpikesV1(Cells)
%
% 
% Example:
%   load Ca2SignalsUdatedAfterViva.mat;
%   VisualiseCellSpikesV1(Cells);

% Set up time axis for movies/rendering
fps = 1/5;
timeaxis = [0:1:143]/fps;

% Get Centroids for ease of processing
 Cr = cat(1,Cells(:).Centroid1);
 Cc = cat(1,Cells(:).Centroid2);
 
 figure(1); % Place the slide cartoon into Figure 1
 plot(Cr,Cc,'o','MarkerFaceColor','b');
 axis ij; % Mode compat with images
 axis off;
 hold on
 for n = 1:length(Cr) % number of cells
     h(n)=render_cells(Cr(n),Cc(n),0,1,1,1,10,[0.8 0.2 0.2]);
     text(Cr(n)-5,Cc(n)+10,20,num2str(n));
 end
 lighting phong;
 hold off;
 
for c = 1:50
    t = SpikeResultsDur(c).locs+1; % t: Array
    MatrixOfPeaksDuring(c,t) = 1; 
end

disp('Press any key to start spike rendering');
for t = 1:100
    for c = 1:50
        if MatrixOfPeaksDuring(c,t)==1
            set(h(c),'FaceColor',[0.8 0.8 0.2]);
        else
            set(h(c),'FaceColor',[0.8 0.2 0.2]);
        end
    end
    pause(0.5);
end
 

 
function h=render_cells(xc,yc,theta,a,b,c,sxyz,Color)
if sxyz ~=0    
    u = 0:0.4:2*pi; % 0.1 gives nice ellipses
    v = (-pi):0.2:(pi);  % 0.1 gives nice ellipses
    x = a*cos(u')*sin(v)*sxyz;
    y = b*sin(u')*sin(v)*sxyz;
    z = c*ones(length(u),1)*cos(v)*sxyz;
    
    %Cm = [ones(length(u),1)]*cos(v+pi/2);
    Cm = [ones(length(u),1)]*cos(v+pi/2);  % experimental
    Cm = 64+round(63*(Cm));  % This ensures we "paint" according to location on ellipse
    %C = Cm;
    C = 1*ones(size(Cm));
    
    x1 =  x*cos(theta) + y*sin(theta);
    y1 = -x*sin(theta) + y*cos(theta);

    K = convhull(x,y,z, 'simplify',true);
    NK = length(K);
   % trisurf(K,x1+xc,y1+yc,z+10,'Facecolor','red','edgecolor','none');
   h=trisurf(K,x1+xc,y1+yc,z+50,'Facecolor',Color,'edgecolor','none');
   %trisurf(K,x1+xc,y1+yc,z+50,C,'edgecolor','none');
   % plot3(x1(K)+xc,y1(K)+yc,z(K)+10,'.')
    axis equal;
    
else
    return
end
 
