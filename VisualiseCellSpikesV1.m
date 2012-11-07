function VisualiseCellSpikesV1(Cells)
% function VisualiseCellSpikesV1(Cells)
%
% 
% Example:
%   load CellsStructure;
%   VisualiseCellSpikesV1(Cells);

% Set up time axis for movies/rendering
fps = 1/5;
timeaxis = [0:1:143]/fps;
nCells = 50;

% Get Centroids for ease of processing
 Cr = cat(1,Cells(:).Centroid1);
 Cc = cat(1,Cells(:).Centroid2);
 
 hf=figure(1); % Place the slide cartoon into Figure 1
 plot(Cr,Cc,'o','MarkerFaceColor','b');
 axis ij; % Mode compat with images
 axis off;
 hold on
 hw = waitbar(0,'Rendering...');
 set(hw,'Color',[1 1 1]); % Black BG
 p = get(hw,'Children');
 hwc=get(p,'Children');
 set(hwc(2),'FaceColor',[1 1 0]);
 set(hwc(2),'EdgeColor',[1 1 0]);
 set(hwc(1),'Color',[1 1 1]);

 for n = 1:length(Cr) % number of cells
     hc(n)=render_cells(Cr(n),Cc(n),0,1,1,1,10,[0.8 0.2 0.2]);
     ht=text(Cr(n)-5,Cc(n)+10,20,num2str(n));
     set(ht,'Color','w');
     waitbar(n/length(Cr),hw)
 end
 delete(hw);
 lighting phong;
 set(hf,'Color',[0,0,0]);
 hold off;
 
MatrixOfPeaksDuring = zeros(nCells,length(timeaxis)); 
for c = 1:nCells
    t = Cells(c).SpikesDuring.time+1; % t: Array
    MatrixOfPeaksDuring(c,t) = 1; 
end

MatrixOfPeaksAfter = zeros(nCells,length(timeaxis)); 
for c = 1:nCells
    t = Cells(c).SpikesAfter.time+1; % t: Array
    MatrixOfPeaksAfter(c,t) = 1; 
end

pause;
disp('Press any key to start spike rendering');
light;

aviobj = avifile('During.avi','FPS',5)
for t = 1:length(timeaxis)
    for c = 1:nCells
        if MatrixOfPeaksDuring(c,t)==1
            set(hc(c),'FaceColor',[0.8 0.8 0.2]);
        else
            set(hc(c),'FaceColor',[0.8 0.2 0.2]);
        end
    end
    SecondString = num2str(timeaxis(t),'%2.1f');
    TimeString = ['t = ',SecondString,'s'];
    htimestamp=text(10,500,[TimeString]);
    set(htimestamp,'Color','w');
    F = getframe;
    aviobj = addframe(aviobj,F);
    pause(0.25);
    delete(htimestamp);
end

aviobj = close(aviobj);
pause;
disp('Press space bar to start spike rendering: After');

for t = 1:length(timeaxis)
    for c = 1:nCells
        if MatrixOfPeaksAfter(c,t)==1
            set(hc(c),'FaceColor',[0.8 0.8 0.2]);
        else
            set(hc(c),'FaceColor',[0.8 0.2 0.2]);
        end
    end
    SecondString = num2str(timeaxis(t),'%2.1f');
    TimeString = ['t = ',SecondString,'s'];
    htimestamp=text(10,510,[TimeString]);
    set(htimestamp,'Color','w');
    pause(0.25);
    delete(htimestamp);
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
 
