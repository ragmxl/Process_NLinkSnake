% 2/21/18
% Modified to accommodate N (4) links...
%
% 1/29/18
% Modified the way the mask is selected (ROI) so that it doesn't include
% the 4 reference markers. Works with the 3 link kinematic snake
%
% *************************************************************************
% ******************************* Mod 1 ***********************************
% 4/19/17
% It works with the 3 Link swimmer with six markers on top, it
% automatically identifies them...
%
% *************************************************************************
% ******************************* Mod 1 ***********************************
% This script allows to automatically process videos of the 3LinkFish, it
% will automatically find 4 markers (on the fish), and detecect the
% appropirate orientation, first maker will be on the head, second on the
% body, third on body and fourth on tail. The way we are able to do this is

% by looking at the first BW frame obtained after processing the color of
% the markers, then we find the pairwise distance between all markers and
% select the four markers that meet the dimensions of the actual 3LinkFish.
% This means that after obtaining the blob properties of the first frame we
% need to apply the transformation to obtain the real coordinates, based on
% the four know reference markers. After the first frame we just carry on
% with the same procedure as before, processing all the frames and it knows
% which markers to track on since we've identified the initial positions.

% *************************************************************************
% ******************************* Mod 2 ***********************************
% This second Modification, will now focus on only processing the area of
% the frame we are interested in. Besides doing the image crop around the
% perimeter of the pool and the mask.
%
% So I need to fully process the first frame and identify the four markers,
% so I can define an area that I have to focus on for the coming frame...

close all
clearvars -except ff        % ff is part of an external for loop
clc  % Comment to use crop and ROI from WorkSpace memory %********%

numberOfLinks = 4;

% rect=round(rect)
% figure(1),clf
%% Load Video, take first frame

% folder='Videos_Spring\';
% ff=5;
% folder='01.13.18\';
folder='4Link\GalaxyS5\';
% folder='4Link\02.23.18\';
file1='20180307_171444';
% file1='4LinkSnake_2_2_Spring_A20_f';
% file1='4LinkSnake_21Spring_A25_f';
% file1='4LinkSnake_f';
% file1='3LinkSnake_3Spring_A20_f';
file2='';
ext='.mp4';
% f=0.1;       % frequency of actuator
f=ff/10             % Use ff in an external for loop
file_middle='';
% if(ff>=10)
%     file_middle=num2str(f*10);
% else if(ff>=1)
%         file_middle=['0',num2str(f*10)];
%     else
%         file_middle=['00',num2str(f*10)];
%     end
% end
file=[folder,file1,file_middle,file2,ext]
% return
video=VideoReader(file);
% video=VideoReader('3Link_LargeEllipse\3LinkSwimmerSymmetric_Spring_f130_A60_S01.mp4');

showTrackingFrames = 0;  % set to "1" to display the region being tracked dynamically, this slows down processing time, "0" to skip
plots=1;    % set to "1" to display plots, "0" for no plots
firstRun=0; % set to "1" if camera or pool/markers have moved, "0" if parameters can be loaded from mat file

if(~firstRun)
    %     load('Preload_Crop_Mask_3LinkFish.mat')
    load('Preload_Crop_Mask_ClassRoom.mat')
        %%% If it hasn't been created use... save('Preload_Crop_Mask_ClassRoom.mat','mask','rect','xrefHigh','yrefHigh')
end

%*** Reference markers ***
% markers=[6;0;6;117.95;232;118.75;232;0];    % measurements in cm (Pool markers)
markers=[0;0;0;173.25;302.2;173.25;302.2;0];        % 3LinkSnake

%*** Select # of markers to track ***
numberOfMarkersOnHighPlane=numberOfLinks*2;

%*** Select color(s) to track ***
%*** HSV Color:
%*** color=1 Red
%*** color=2 Green
%*** color=3 Blue
%*** color=4 Yellow

color=4; % could use more colors ex: color=[1,2,3,4]
startFrame=1;
im=read(video,startFrame);


%% Crop area of interest.
if(firstRun)
    display('Select area to crop')
    [im,rect]=imcrop(im);             % Select crop area        %********%
else
    im=imcrop(im,rect);                 % Crop area in memory %********%
end
im1cropped=im;

%% Create Polygon for ROI to create mask, or load one from file.
if(firstRun)
    display('Create polygon for ROI')
    mask=roipoly(im);              % Select ROI or use from memory %********%
end

%% Process Frames
P=5;    % points used for bwareaopen(BW,P,CONN) % Used 8 for in-lab videos
CONN=4; % connectivity

P_ref=4;        % Same for reference markers
CONN_ref=4;

t0=tic;
endFrame=100;
fps=video.frame;
% N=video.NumberOfFrames-(startFrame-1);   %This or ... N = last frame I want
N=660;
N=400;
s(N)=struct('centroids',[]);
bHigh=zeros(N,2,numberOfMarkersOnHighPlane);
for i=1:N,    %i:N
    if(i>1)
        im=imcrop(read(video,i+(startFrame-1)),rect2);
        BW=detectColorHSV(im,color);  
        P=1;
    else
        % _refMarkers to process at least frame 15 for the fixed markers,
        % sometimes if we process the initial frame, there is some issue
        % detecting the color, the image is still adjusting...
        im_refMarkers=imcrop(read(video,startFrame+15),rect);
        % BW_refMarkers=detectColorHSV(im_refMarkers.*uint8(repmat(mask,[1,1,3])),color);
        % With the previous line commented I am not using a mask for the 4 ref markers...
        BW_refMarkers=detectColorHSV(im_refMarkers,color);
        BW_refMarkers=bwareaopen(BW_refMarkers,P_ref,CONN_ref);
        % We also do it without _refMarkers, for normal processing of
        % markers on robot...
        BW=detectColorHSV(im.*uint8(repmat(mask,[1,1,3])),color);
        BW=bwareaopen(BW,P,CONN);
        if(firstRun)
            display('Select the 4 reference points for the High Plane,')
            figure(1),clf
            imshow(im)
            [xrefHigh,yrefHigh]=ginput(4);    % comment to use previous points (don't clear)
        end
        if(plots)
            figure(1),clf
            imshow(im)
            title('First frame');
            figure(2),clf
            imshow(BW,'InitialMagnification','fit')
            title('BW first frame');
            figure(3),clf
            imshow(BW_refMarkers,'InitialMagnification','fit')
            title('BW first frame for reference markers, it really is frame 16');
        end
    end
   
    % Blob Properties
    if(i~=1)
        BW=bwareaopen(BW,P,CONN);     % Remove regions with a small area...
    end
    
    % % To Test orientation identification of the snake, we rotate image 180
    % % degrees
    % figure(10),clf
    % imshow(BW)
    % figure(11),clf
    % BW=imrotate(BW,180);
    % imshow(BW)
    
    blob=regionprops(BW,'centroid');
    centroids=cat(1,blob.Centroid);
    
    if(i>1)
        % Give the centroids in the coordinates of the whole experiment area
        centroids(:,1)=centroids(:,1)+xOffset-1;
        centroids(:,2)=centroids(:,2)+yOffset-1;
        % Track the closest marker
        for k=1:numberOfMarkersOnHighPlane,
            bHigh(i,:,k)=closestBall(centroids,bHigh(i-1,:,k));
        end

        % Define bounds of region to track...
        oldMinX=min(bHigh(i,1,:));
        oldMaxX=max(bHigh(i,1,:));
        oldMinY=min(bHigh(i,2,:));
        oldMaxY=max(bHigh(i,2,:));
        xOffset=round(oldMinX)-50;              % May need to adapt size larger for the Large Ellipses
        yOffset=round(oldMinY)-50;              % was 50
        xLength=round(oldMaxX-oldMinX)+100;     % and 100
        yLength=round(oldMaxY-oldMinY)+100;
        rect2=[rect(1)-1+xOffset,rect(2)-1+yOffset,xLength,yLength];
    end
    
    s(i).centroids=centroids;
    if(showTrackingFrames)
        if(mod(i,10)==0)
            display(['Processed Frames: ',num2str(i)])
            figure(23)
            imshow(imcrop(read(video,i),rect2))
            pause(0.1)
        end
    end
    
    if(i==1)
        % _refMarkers to process at least frame 15 for the fixed markers,
        % sometimes if we process the initial frame, there is some issue
        % detecting the color, the image is still adjusting...
        blob_refMarkers=regionprops(BW_refMarkers,'centroid');
        centroids_refMarkers=cat(1,blob_refMarkers.Centroid);
        % we still do blob and centroids without _refMarkers, to process
        % first frame for the markers on the robot...
        blob=regionprops(BW,'centroid');
        centroids=cat(1,blob.Centroid);
        ref1High=closestBall(centroids_refMarkers,[xrefHigh(1),yrefHigh(1)]);
        ref2High=closestBall(centroids_refMarkers,[xrefHigh(2),yrefHigh(2)]);
        ref3High=closestBall(centroids_refMarkers,[xrefHigh(3),yrefHigh(3)]);
        ref4High=closestBall(centroids_refMarkers,[xrefHigh(4),yrefHigh(4)]);
        
        % Move origin to lower left corner
        %******** Dimensions of cropped image
        dim=size(im1cropped);
        xmax=dim(2);
        ymax=dim(1);
        
        %******** Flip coordinates to have the origin in the lower left corner
        s1Flipped=s(1).centroids;
        s1Flipped(:,2)=-s1Flipped(:,2)+ymax;
        ref1High(2)=-ref1High(2)+ymax;
        ref2High(2)=-ref2High(2)+ymax;
        ref3High(2)=-ref3High(2)+ymax;
        ref4High(2)=-ref4High(2)+ymax;
        
        % Projective Transformation (camera coordinates to inertial frame)
        %Find Transformation
        %*** High Plane ***   (1. Bottom Left, 2. Top Left, 3. Top Right, 4. Bottom Right)
        THigh=createProjectiveTransformation([ref1High(1);ref1High(2);ref2High(1);ref2High(2);...
            ref3High(1);ref3High(2);ref4High(1);ref4High(2)],...
            1*markers);
        
        %Projection
        %*** High Plane ***
        s1Real=projectCoordinates(s1Flipped,THigh);
        
        %Now we calculate the pairwise distances
        nBlobs=length(s1Real);
        dist=zeros(nBlobs);
        for k=1:nBlobs-1,
            for j=k+1:nBlobs,
                dist(k,j)=norm(s1Real(k,:)-s1Real(j,:));
            end
        end,    dist
        
        %HEAD
        headDist=8.7;
        headMatrix=abs(dist-headDist);
        [head_i,head_j]=find(headMatrix==min(headMatrix(:)));
        head_centroid=(s1Real(head_i,:)+s1Real(head_j,:))/2;
        %BODY LINKS
        nBodyLinks=numberOfLinks-1;       % For the 4 Link Snake this should be 3
        bodyLinkIJ=zeros(nBodyLinks,2); % this will indicate the two blobs per link, blob I and blob J
        bodyDist=6.5;
        bodyMatrix=abs(dist-bodyDist);
        for k=1:nBodyLinks,
            [bodyLinkIJ(k,1),bodyLinkIJ(k,2)]=find(bodyMatrix==min(bodyMatrix(:)));
            bodyMatrix(bodyLinkIJ(k,1),bodyLinkIJ(k,2))=k*1111;
        end
        
        % So I know which pairs represent the Head and Body Links
        % I should now order them from Head to tail...
        
        head_i_body_dist=norm(s1Real(head_i,:)-s1Real(bodyLinkIJ(1,1),:));     % obtain distance to any blob from any other link, they should start in a straight line...
        head_j_body_dist=norm(s1Real(head_j,:)-s1Real(bodyLinkIJ(1,1),:));     % obtain distance to any blob from any other link, they should start in a straight line...
        if( head_i_body_dist > head_j_body_dist )
            head1_Blob = head_i;
            head2_Blob = head_j;
        else
            head1_Blob = head_j;
            head2_Blob = head_i;
        end
        
        % I know now the order of the blobs for the Head link
        % now I need to order the rest of body links...
        % for body1, body2, body3 ...
        % I will first creat a matrix with distances from all blobs
        % of interest to the HEAD Centroid and order therm from
        % closest to farthest...
        
        bodyBlobs=zeros(2*nBodyLinks,1);        % Put all blobls indices in a single vector
        for k=1:nBodyLinks,
            bodyBlobs(2*(k-1)+1:2*k)=bodyLinkIJ(k,:)';
        end
        
        bodyBlobsDistances=zeros(2*nBodyLinks,1);
        for k=1:2*nBodyLinks,
            bodyBlobsDistances(k)=norm(s1Real(bodyBlobs(k))-s1Real(head2_Blob));
        end
        
        bodyBlobsAndDistances=sortrows([bodyBlobsDistances,bodyBlobs]);
        bodyBlobsOrdered=[head1_Blob;head2_Blob;bodyBlobsAndDistances(:,2)];
        
        % Now we know the order of all markers...
        % (bodyBlobsOrdered)... this includes the head blobs
        fprintf('\n\tOrdered blobs, Head to Tail: ')
        fprintf('%d ',[bodyBlobsOrdered'])
        fprintf('\n\n')
        
        %Now we already now our starting markers positions... 2N markers
        %for an N Link fish/snake... (this will be the initial markers)
        for k=1:numberOfMarkersOnHighPlane,
            bHigh(1,:,k)=s(1).centroids(bodyBlobsOrdered(k),:);
        end
        
        oldMinX=min(bHigh(1,1,:));
        oldMaxX=max(bHigh(1,1,:));
        oldMinY=min(bHigh(1,2,:));
        oldMaxY=max(bHigh(1,2,:));
        xOffset=round(oldMinX)-50;
        yOffset=round(oldMinY)-50;
        xLength=round(oldMaxX-oldMinX)+100;
        yLength=round(oldMaxY-oldMinY)+100;
        rect2=[rect(1)-1+xOffset,rect(2)-1+yOffset,xLength,yLength];
        %         return
    end
    
end     % bHigh 1 through 6 contains the positions of the six markers of the swimmer
display(['Total Frames: ',num2str(i)])
display(['Total Processing time: ',num2str(toc(t0)),' seconds.'])

imLastCropped=imcrop(read(video,N),rect);

%%
if(plots)
    figure(1)
    hold on
    colorstring='brgmcykbrgmcykbrgmcyk';
    for i=1:numberOfMarkersOnHighPlane
        plot(bHigh(:,1,i),bHigh(:,2,i),'Color',colorstring(i))
    end
    figure(4),clf
    imshow(imLastCropped,'InitialMagnification','fit')
    title('Last frame');
    hold on
    for i=1:numberOfMarkersOnHighPlane
        plot(bHigh(:,1,i),bHigh(:,2,i),'Color',colorstring(i))
    end
end

%% Move origin to lower left corner
%******** Dimensions of cropped image
dim=size(im1cropped);
xmax=dim(2);
ymax=dim(1);

%******** Flip coordinates to have the origin in the lower left corner
for i=1:numberOfMarkersOnHighPlane
    bbHigh(:,:,i)=bHigh(:,:,i);
end

bbHigh(:,2,:)=-bbHigh(:,2,:)+ymax;


%% Projective Transformation (camera coordinates to inertial frame)
%Projection
%*** High Plane ***
bHighReal=zeros(N,2,numberOfMarkersOnHighPlane);

for i=1:numberOfMarkersOnHighPlane
    bHighReal(:,:,i)=projectCoordinates(bbHigh(:,:,i),THigh);
end
refReal=projectCoordinates([ref1High;ref2High;ref3High;ref4High;ref1High],THigh)

if(plots)
    figure(5),clf
    hold on
    for i=1:numberOfMarkersOnHighPlane
        plot(bHighReal(:,1,i),bHighReal(:,2,i),'Color',colorstring(i))
    end
    plot((bHighReal(:,1,3)+bHighReal(:,1,4))/2,(bHighReal(:,2,3)+bHighReal(:,2,4))/2,'Color',colorstring(i+1))
    xlabel('x (cm)')
    ylabel('y (cm)')
    title('Real coordinates');
    % for i=1:numberOfMarkersOnLowPlane
    %     plot(bLowReal(:,1,i),bLowReal(:,2,i),'Color',colorstring(numberOfMarkersOnHighPlane+1+i))
    % end
    plot(refReal(:,1),refReal(:,2))
    axis equal tight
end

%% Plot fixed length between markers
if(plots)
    figure(6),clf
    for k=1:numberOfLinks,
        diff=bHighReal(:,:,2*(k-1)+1)-bHighReal(:,:,2*k);
        for i=1:length(diff),
            d(i)=norm(diff(i,:));
        end
        subplot(numberOfLinks,1,k)
        plot(d)
        hold on
        if(k==1)
            plot([0,length(diff)],[headDist,headDist],'k');
        else
            plot([0,length(diff)],[bodyDist,bodyDist],'k');
        end
        plot([0,length(diff)],[mean(d),mean(d)],':k');
        legend('raw','real','mean')
        title(['Link ',num2str(k),' markers'' fixed length'])
        %     xlabel('frame')
        ylabel('Length (cm)')
    end
    xlabel('frame')
end

% return

%% Save data...
%save('savedData\3LinkSwimmerSymmetric_Spring_f140_A60_S01.mat','fps','bHighReal','f','im1cropped','imLastCropped','bHigh')
save([folder,'savedData\',file1,file_middle,file2,'.mat'],'fps','bHighReal','f','im1cropped','imLastCropped','bHigh')
% save(['savedData\',file1,file_middle,file2,'.mat'],'fps','bHighReal','f','im1cropped','imLastCropped','bHigh')
% PostProcess_3LinkEllipse
figure(6)   % See errors...
set (figure(6), 'Units', 'normalized', 'Position', [0,0,1,1]);
% load handel
% sound(y,Fs)
figure(6), pause(1),
% saveAsPdf(6,['Processed_MarkersDistances\FixedMarkers_',file1,file_middle,file2])
saveAsPdf(6,[folder,'Processed_MarkersDistances\FixedMarkers_',file1,file_middle,file2])
pause(3)
