function VisualiseCellSpikesV1(Cells,makeavi)
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

if ~exist('makeavi','var')
    makeavi = 0;
end

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

if makeavi
    aviobj = avifile('DuringAfter.avi','FPS',5)
end
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
    htWhen = text(10,500,'During Flow');
    htimestamp=text(10,530,[TimeString]);
    set(htWhen,'Color','w');
    set(htimestamp,'Color','w');
    
    
    %%
    if makeavi
        hf = gcf;
        F = getframe(hf,[100,1,500,500]);
        aviobj = addframe(aviobj,F);
    end
    %% For debugging only...
    % figure(2);image(F.cdata);
    % pause
    pause(0.25);% If viewing in RT, remove "%"
    delete(htWhen);
    delete(htimestamp);
end


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
    
    htWhen = text(10,500,'Flow Cessation');
    htimestamp=text(10,530,[TimeString]);
    set(htWhen,'Color','w');
    set(htimestamp,'Color','w');
    
    if makeavi
        hf = gcf;
        F = getframe(hf,[100,1,500,500]);
        aviobj = addframe(aviobj,F);
    end
    
    %pause(0.25); % If viewing in RT, remove "%"
    delete(htWhen);
    delete(htimestamp);
end

if makeavi
    aviobj = close(aviobj);
end

 
function h=render_cells(xc,yc,theta,a,b,c,sxyz,Color)
if sxyz ~=0    
    u = 0:0.4:2*pi; % 0.1 gives smoother ellipses
    v = (-pi):0.2:(pi);  % 0.1 gives smoother result
    x = a*cos(u')*sin(v)*sxyz;
    y = b*sin(u')*sin(v)*sxyz;
    z = c*ones(length(u),1)*cos(v)*sxyz;
    
    x1 =  x*cos(theta) + y*sin(theta);
    y1 = -x*sin(theta) + y*cos(theta);

    K = convhull(x,y,z, 'simplify',true);
    NK = length(K);
   
    h=trisurf(K,x1+xc,y1+yc,z+50,'Facecolor',Color,'edgecolor','none');
   
    axis equal;
    
else
    return
end
 
