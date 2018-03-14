% 2/21/18
% Modified to accommodate 4 links...
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
% rect=round(rect)
% figure(1),clf
%% Load Video, take first frame

% folder='Videos_Spring\';
folder='01.13.18\';
file1='3LinkSnake_3Spring_A20_f';
% file1='3LinkSwimmer_Symmetric_FullyActuated_f';
% file1='3LinkSwimmer_Symmetric_Spring_f';
% file1='3LinkSwimmer_FatHead_Spring_f';
% file1='3LinkSwimmer_FatHead_EqSpacing_Spring_f';
% file1='3LinkSwimmer_FatHead_EqLeverage_Spring_f';
% file1='3LinkSwimmer_SkinnyHead_Spring_f';
% file1='3LinkSwimmer_SkinnyHead_EqSpacing_Spring_f';
% file1='3LinkSwimmer_SkinnyHead_EqLeverage_Spring_f';
% file1='testAngles3';
% file2='_A60_S01';
% file2='_A30_S02';
% file2='_Test05';
% file2='_A60_S03';
file2='';
ext='.mp4';
% f=0.1;       % frequency of actuator
f=ff/10             % Use ff in an external for loop
file_middle='';
if(ff>=10)
    file_middle=num2str(f*10);
else if(ff>=1)
        file_middle=['0',num2str(f*10)];
    else
        file_middle=['00',num2str(f*10)];
    end
end
file=[folder,file1,file_middle,file2,ext]
% return
video=VideoReader(file);
% video=VideoReader('3Link_LargeEllipse\3LinkSwimmerSymmetric_Spring_f130_A60_S01.mp4');

plots=1;    % set to "1" to display plots, "0" for no plots
firstRun=1; % set to "1" if camera or pool/markers have moved, "0" if parameters can be loaded from mat file

if(~firstRun)
%     load('Preload_Crop_Mask_3LinkFish.mat')
    load('Preload_Crop_Mask_3LinkSnake.mat')
end

%*** Reference markers ***
% markers=[6;0;6;117.95;232;118.75;232;0];    % measurements in cm
markers=[0;0;0;135.3;231.7;135.3;230.6;0];        % 3LinkSnake

%*** Select # of markers to track ***
numberOfMarkersOnHighPlane=6;
numberOfMarkersOnLowPlane=0;

%*** Select color(s) to track ***
%*** HSV Color:
%*** color=1 Red
%*** color=2 Green
%*** color=3 Blue
%*** color=4 Yellow
% color=4;
color=4; % could use more colors ex: color=[1,2,3,4]

% t1a=tic;
% startFrame=28;    %hydrox-floor-trimmed2.mov
% startFrame=69;      %20150807_140824.mp4
% startFrame=95;      %2balls_Red_and_Blue.mp4
% startFrame=75;
startFrame=1;
im=read(video,startFrame);
% im=readFrame(video);
% t1b=toc(t1a);
% display(num2str(t1b))
% figure(1)
% imshow(im,'InitialMagnification','fit')

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

%% ************* Mod 2 *************
% Get transformation ready, from the 4 known reference markers
% % Select 4 Reference points
% %*** High Plane ***
% display('Select the 4 reference points for the High Plane,')
% figure(4),clf
% % im=imcrop(read(video,200),rect);
% im=imcrop(read(video,1),rect);
% imshow(im);
% % [xrefHigh,yrefHigh]=ginput(4);    % comment to use previous points (don't clear)
% % display([x,y])
% BW=detectColorHSV(im.*uint8(repmat(mask,[1,1,3])),4); % Yellow color
% clf
% imshow(BW)

%% Process Frames
t0=tic;
endFrame=100;
fps=video.frame;
N=video.NumberOfFrames-(startFrame-1);   %This or ...
% N=100;                                          %This
% N=endFrame-(startFrame-1);
% N=100;
s(N)=struct('centroids',[]);
% for i=1:video.NumberOfFrames,
bHigh=zeros(N,2,numberOfMarkersOnHighPlane);
for i=1:N,    %i:N
    if(i>1)
% %         im=imcrop(read(video,i+(startFrame-1)),rect);
        im=imcrop(read(video,i+(startFrame-1)),rect2);
%         im=imcrop(readFrame(video),rect2);
        BW=detectColorHSV(im,color);  %
    else
        % _refMarkers to process at least frame 15 for the fixed markers,
        % sometimes if we process the initial frame, there is some issue
        % detecting the color, the image is still adjusting...
        im_refMarkers=imcrop(read(video,startFrame+15),rect);
        % BW_refMarkers=detectColorHSV(im_refMarkers.*uint8(repmat(mask,[1,1,3])),color);
        % With the previous line commented I am not using a mask for the 4 ref markers...
        BW_refMarkers=detectColorHSV(im_refMarkers,color);
        BW_refMarkers=bwareaopen(BW_refMarkers,8,4);
        % We also do it without _refMarkers, for normal processing of
        % markers on robot...
        BW=detectColorHSV(im.*uint8(repmat(mask,[1,1,3])),color);
        BW=bwareaopen(BW,8,4);
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

% if(i==1)
%     figure(2)%5
%     BW=bwareaopen(BW,8,4);    
%     imshow(BW,'InitialMagnification','fit')
% end
%********************************************************************%
%********************************************************************%

% Blob Properties
if(i~=1)
    BW=bwareaopen(BW,8,4);     %Used to be (BW,18) (BW,6,4)      % Remove regions with a small area...
end
blob=regionprops(BW,'centroid');
centroids=cat(1,blob.Centroid);

if(i>1)
    % Give the centroids in the coordinates of the whole experiment area
    centroids(:,1)=centroids(:,1)+xOffset-1;
    centroids(:,2)=centroids(:,2)+yOffset-1;
    % Track the closest marker
    bHigh(i,:,1)=closestBall(centroids,bHigh(i-1,:,1));
    bHigh(i,:,2)=closestBall(centroids,bHigh(i-1,:,2));
    bHigh(i,:,3)=closestBall(centroids,bHigh(i-1,:,3));
    bHigh(i,:,4)=closestBall(centroids,bHigh(i-1,:,4));
    bHigh(i,:,5)=closestBall(centroids,bHigh(i-1,:,5));
    bHigh(i,:,6)=closestBall(centroids,bHigh(i-1,:,6));
    %
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
%     if(mod(i,10)==0)
%         display(['Processed Frames: ',num2str(i)])
% %         figure(23)
% %         imshow(imcrop(read(video,i),rect2))
%     end
    
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
%         figure(10),clf
%         imshow(im)
%         figure(11),clf
%         imshow(bwareaopen(BW,8,4))
        
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
        end
        
        %Now we identify the six markers on the 3-Link Swimmer
%         distance of tail markers:  7.62 cm
%         distance of body markers: 19.05 cm
%         distance of headmarkers: 11.43 cm 

%         tailDist=7.62;
        tailDist = 10.3;
        tailMatrix=abs(dist-tailDist);
        [tail_i,tail_j]=find(tailMatrix==min(tailMatrix(:)));
        tail_centroid=(s1Real(tail_i,:)+s1Real(tail_j,:))/2;

%         bodyDist=19.05;
        bodyDist=6.5;
        bodyMatrix=abs(dist-bodyDist);
        [body_i,body_j]=find(bodyMatrix==min(bodyMatrix(:)));
        body_centroid=(s1Real(body_i,:)+s1Real(body_j,:))/2;
        
%         headDist=11.43;
        headDist=8.7;
        headMatrix=abs(dist-headDist);
        [head_i,head_j]=find(headMatrix==min(headMatrix(:)));
        head_centroid=(s1Real(head_i,:)+s1Real(head_j,:))/2;
        
        % So I know which pairs represent the Head, Body and Tail
        % I should now order them from Head to tail...
        
        head_i_body_dist=norm(s1Real(head_i,:)-body_centroid);
        head_j_body_dist=norm(s1Real(head_j,:)-body_centroid);
        if( head_i_body_dist > head_j_body_dist )
            head1_Blob = head_i;
            head2_Blob = head_j;
        else
            head1_Blob = head_j;
            head2_Blob = head_i;
        end
        
        body_i_head_dist=norm(s1Real(body_i,:)-head_centroid);
        body_j_head_dist=norm(s1Real(body_j,:)-head_centroid);
        if( body_i_head_dist < body_j_head_dist )
            body3_Blob = body_i;
            body4_Blob = body_j;
        else
            body3_Blob = body_j;
            body4_Blob = body_i;
        end
        
        tail_i_body_dist=norm(s1Real(tail_i,:)-body_centroid);
        tail_j_body_dist=norm(s1Real(tail_j,:)-body_centroid);
        if( tail_i_body_dist < tail_j_body_dist )
            tail5_Blob = tail_i;
            tail6_Blob = tail_j;
        else
            tail5_Blob = tail_j;
            tail6_Blob = tail_i;
        end
        
        % So now I know the order of the 6 blobs.
                                                                                            % ***************************************
                                                                                            % This was for the other geometry markers
%         body_and_head_Dist=11.43;
%         bodyMatrix=abs(dist-body_and_head_Dist);
%         [temp1_i,temp1_j]=find(bodyMatrix==min(bodyMatrix(:)));
%         bodyMatrix(temp1_i,temp1_j)=10000;
%         [temp2_i,temp2_j]=find(bodyMatrix==min(bodyMatrix(:)));
%         temp1_centroid=(s1Real(temp1_i,:)+s1Real(temp1_j,:))/2;
%         temp2_centroid=(s1Real(temp2_i,:)+s1Real(temp2_j,:))/2;
%         
%         dist_temp1=norm(temp1_centroid-tail_centroid);
%         dist_temp2=norm(temp2_centroid-tail_centroid);
%         
%         if(dist_temp1 > dist_temp2)     %means temp 2 is the body
%             head1_temp=s1Real(temp1_i,:);
%             head2_temp=s1Real(temp1_j,:);
%             head_i=temp1_i;
%             head_j=temp1_j;
%             body1_temp=s1Real(temp2_i,:);
%             body2_temp=s1Real(temp2_j,:);
%             body_i=temp2_i;
%             body_j=temp2_j;
%         else
%             head1_temp=s1Real(temp2_i,:);
%             head2_temp=s1Real(temp2_j,:);
%             head_i=temp2_i;
%             head_j=temp2_j;
%             body1_temp=s1Real(temp1_i,:);
%             body2_temp=s1Real(temp1_j,:);
%             body_i=temp1_i;
%             body_j=temp1_j;
%         end        
%         head1_tail_dist=norm(head1_temp-tail_centroid);
%         head2_tail_dist=norm(head2_temp-tail_centroid);
%         body1_tail_dist=norm(body1_temp-tail_centroid);
%         body2_tail_dist=norm(body2_temp-tail_centroid);
%             
%         if(head1_tail_dist > head2_tail_dist)   %means 1 is the first blob
% %             head1_Blob=head1_temp;
% %             head2_Blob=head2_temp;
%             head1_Blob=head_i;
%             head2_Blob=head_j;
%         else
% %             head1_Blob=head2_temp;
% %             head2_Blob=head1_temp;
%             head1_Blob=head_j;
%             head2_Blob=head_i;
%         end
%         if(body1_tail_dist > body2_tail_dist)   %means 1 is the third blob
% %             body3_Blob=body1_temp;
% %             body4_Blob=body2_temp;
%             body3_Blob=body_i;
%             body4_Blob=body_j;
%         else
% %             body3_Blob=body2_temp;
% %             body4_Blob=body1_temp;
%             body3_Blob=body_j;
%             body4_Blob=body_i;
%         end
%         %Now find which tail 5 and 6
%         body_centroid=(s1Real(body_i,:)+s1Real(body_j,:))/2;
%         tail1_body_dist=norm(s1Real(tail_i,:)-body_centroid);
%         tail2_body_dist=norm(s1Real(tail_j,:)-body_centroid);
%         if(tail1_body_dist > tail2_body_dist)   %means j is the fifth blob
%             tail5_Blob=tail_j;
%             tail6_Blob=tail_i;
%         else
%             tail5_Blob=tail_i;
%             tail6_Blob=tail_j;
%         end               
                                                                                            % ***************************************
                                                                                            % This was for the other geometry markers

%         return
        
%         headBodyDist=8.8011;    
%         headMatrix=abs(dist-headBodyDist);
%         [hb_i,hb_j]=find(headMatrix==min(headMatrix(:)));
%         bodyBodyDist=7.874;
%         bodyMatrix=abs(dist-bodyBodyDist);
%         [bb_i,bb_j]=find(bodyMatrix==min(bodyMatrix(:)));
%         tailBodyDist=6.8453;
%         tailMatrix=abs(dist-tailBodyDist);
%         [tb_i,tb_j]=find(tailMatrix==min(tailMatrix(:)));
        
        %We now know which blobs belong to the fish, now we need to find out which
%         %is which
%         if(hb_i == bb_i || hb_i == bb_j)
%             headBlob = hb_j;
%             bodyHeadBlob = hb_i;
%         else
%             headBlob = hb_i;
%             bodyHeadBlob = hb_j;
%         end
%         
%         if(tb_i == bb_i || tb_i == bb_j)
%             tailBlob = tb_j;
%             bodyTailBlob = tb_i;
%         else
%             tailBlob = tb_i;
%             bodyTailBlob = tb_j;
%         end
%         
%         return
        
        %Now we already now our starting markers positions... 4 markers for fish
        bHigh(1,:,1)=s(1).centroids(head1_Blob,:);
        bHigh(1,:,2)=s(1).centroids(head2_Blob,:);
        bHigh(1,:,3)=s(1).centroids(body3_Blob,:);
        bHigh(1,:,4)=s(1).centroids(body4_Blob,:);
        bHigh(1,:,5)=s(1).centroids(tail5_Blob,:);
        bHigh(1,:,6)=s(1).centroids(tail6_Blob,:);
        oldMinX=min(bHigh(1,1,:));
        oldMaxX=max(bHigh(1,1,:));
        oldMinY=min(bHigh(1,2,:));
        oldMaxY=max(bHigh(1,2,:));
        xOffset=round(oldMinX)-50;
        yOffset=round(oldMinY)-50;
        xLength=round(oldMaxX-oldMinX)+100;
        yLength=round(oldMaxY-oldMinY)+100;
        rect2=[rect(1)-1+xOffset,rect(2)-1+yOffset,xLength,yLength];
    end
        
    
end     % bHigh 1 through 6 contains the positions of the six markers of the swimmer
display(['Total Frames: ',num2str(i)])
display(['Total Processing time: ',num2str(toc(t0)),' seconds.'])

imLastCropped=imcrop(read(video,N),rect);
% return
%% Track balls (markers)
% % % % % % ***************** Modification *************************
% % % % % % ********************************************************
% % % % % % To automatically detect the four makers that make up the fish, detecting
% % % % % % head, body and tail. I will need to first get the real coordinates of
% % % % % % every blob, and then compare the pairwise distances to see which pairs of
% % % % % % blobs belong to the body. The distances of the markers are known so I can
% % % % % % detect them that way.
% % % % % 
% % % % % % Select 4 Reference points
% % % % % %*** High Plane ***
% % % % % display('Select the 4 reference points for the High Plane,')
% % % % % figure(4),clf
% % % % % % im=imcrop(read(video,200),rect);
% % % % % im=imcrop(read(video,1),rect);
% % % % % imshow(im);
% % % % % % [xrefHigh,yrefHigh]=ginput(4);    % comment to use previous points (don't clear)
% % % % % % display([x,y])
% % % % % % % BW=detectColorHSV(im.*uint8(repmat(mask,[1,1,3])),4); % Yellow color
% % % % % % % clf
% % % % % % % imshow(BW)
% % % % % % % blob=regionprops(BW,'centroid');
% % % % % % % centroids=cat(1,blob.Centroid);
% % % % % % % ref1High=closestBall(centroids,[xrefHigh(1),yrefHigh(1)]);
% % % % % % % ref2High=closestBall(centroids,[xrefHigh(2),yrefHigh(2)]);
% % % % % % % ref3High=closestBall(centroids,[xrefHigh(3),yrefHigh(3)]);
% % % % % % % ref4High=closestBall(centroids,[xrefHigh(4),yrefHigh(4)]);
% % % % % % % figure(10),clf
% % % % % % % imshow(im)
% % % % % % % figure(11),clf
% % % % % % % imshow(bwareaopen(BW,8,4))
% % % % % 
% % % % % % Move origin to lower left corner
% % % % % %******** Dimensions of cropped image
% % % % % dim=size(im1cropped);
% % % % % xmax=dim(2);
% % % % % ymax=dim(1);
% % % % % 
% % % % % %******** Flip coordinates to have the origin in the lower left corner
% % % % % s1Flipped=s(1).centroids;
% % % % % s1Flipped(:,2)=-s1Flipped(:,2)+ymax;
% % % % % ref1High(2)=-ref1High(2)+ymax;
% % % % % ref2High(2)=-ref2High(2)+ymax;
% % % % % ref3High(2)=-ref3High(2)+ymax;
% % % % % ref4High(2)=-ref4High(2)+ymax;
% % % % % 
% % % % % % Projective Transformation (camera coordinates to inertial frame)
% % % % % %Find Transformation
% % % % % %*** High Plane ***   (1. Bottom Left, 2. Top Left, 3. Top Right, 4. Bottom Right)
% % % % % % % THigh=createProjectiveTransformation([ref1High(1);ref1High(2);ref2High(1);ref2High(2);...
% % % % % % %     ref3High(1);ref3High(2);ref4High(1);ref4High(2)],...
% % % % % % %     1*[6;0;6;117.95;232;118.75;232;0]);
% % % % % % %already obtained when processed first frame...
% % % % % 
% % % % % %Projection
% % % % % %*** High Plane ***
% % % % % s1Real=s1Flipped;
% % % % % s1Real=projectCoordinates(s1Flipped,THigh);
% % % % % 
% % % % % %Now we calculate the pairwise distances
% % % % % % nBlobs=length(s1Real);
% % % % % % dist=zeros(nBlobs);
% % % % % % for i=1:nBlobs-1,
% % % % % %     for j=i+1:nBlobs,
% % % % % %         dist(i,j)=norm(s1Real(i,:)-s1Real(j,:));
% % % % % %     end
% % % % % % end
% % % % % % headBodyDist=8.8011;
% % % % % % headMatrix=abs(dist-headBodyDist);
% % % % % % [hb_i,hb_j]=find(headMatrix==min(headMatrix(:)))
% % % % % % bodyBodyDist=7.874;
% % % % % % bodyMatrix=abs(dist-bodyBodyDist);
% % % % % % [bb_i,bb_j]=find(bodyMatrix==min(bodyMatrix(:)))
% % % % % % tailBodyDist=6.8453;
% % % % % % tailMatrix=abs(dist-tailBodyDist);
% % % % % % [tb_i,tb_j]=find(tailMatrix==min(tailMatrix(:)))
% % % % % % 
% % % % % % %We now know which blobs belong to the fish, now we need to find out which
% % % % % % %is which
% % % % % % if(hb_i == bb_i || hb_i == bb_j)
% % % % % %     headBlob = hb_j;
% % % % % %     bodyHeadBlob = hb_i;
% % % % % % else
% % % % % %     headBlob = hb_i;
% % % % % %     bodyHeadBlob = hb_j;
% % % % % % end
% % % % % % 
% % % % % % if(tb_i == bb_i || tb_i == bb_j)
% % % % % %     tailBlob = tb_j;
% % % % % %     bodyTailBlob = tb_i;
% % % % % % else
% % % % % %     tailBlob = tb_i;
% % % % % %     bodyTailBlob = tb_j;
% % % % % % end
% % % % % % 
% % % % % % %Now we already now our starting markers positions... 4 markers for fish
% % % % % % bHigh=zeros(N,2,numberOfMarkersOnHighPlane);
% % % % % % bHigh(1,:,1)=s(1).centroids(headBlob,:);
% % % % % % bHigh(1,:,2)=s(1).centroids(bodyHeadBlob,:);
% % % % % % bHigh(1,:,3)=s(1).centroids(bodyTailBlob,:);
% % % % % % bHigh(1,:,4)=s(1).centroids(tailBlob,:);
% % % % % 
% % % % % % ***************** Modification *************************
% % % % % % ********************************************************

%%
%*** High Plane ***
%******Commented to use modification of auto-detect 4 markers *************
% display(['Select the initial ',num2str(numberOfMarkersOnHighPlane),' markers to track on the High Plane,'])
% figure(1),clf
% imshow(im1cropped,'InitialMagnification','fit')
% [xtr,ytr]=ginput(numberOfMarkersOnHighPlane);    % comment to use previous points (don't clear)
% display([xtr,ytr])
% 
% bHigh=zeros(N,2,numberOfMarkersOnHighPlane);
% for i=1:numberOfMarkersOnHighPlane
%     bHigh(1,:,i)=closestBall(s(1).centroids,[xtr(i),ytr(i)]);
% end
%******Commented to use modification of auto-detect 4 markers *************
% % for i=2:N
% %     for j=1:numberOfMarkersOnHighPlane
% %         bHigh(i,:,j)=closestBall(s(i).centroids,bHigh(i-1,:,j));
% %     end
% % end
%May be redundant, already done when processed frames...

if(plots)
    figure(1)
    hold on
    colorstring='brgmcykbrg';
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
% pause
%*** Low Plane ***
% % display(['Select the initial ',num2str(numberOfMarkersOnLowPlane),' marker(s) to track on the Low Plane,'])
% % figure(1),clf
% % imshow(im1cropped,'InitialMagnification','fit')
% % [xtr,ytr]=ginput(numberOfMarkersOnLowPlane);    % comment to use previous points (don't clear)
% % display([xtr,ytr])
% % 
% % bLow=zeros(N,2,numberOfMarkersOnLowPlane);
% % for i=1:numberOfMarkersOnLowPlane
% %     bLow(1,:,i)=closestBall(s(1).centroids,[xtr(i),ytr(i)]);
% % end
% % for i=2:N
% %     for j=1:numberOfMarkersOnLowPlane
% %         bLow(i,:,j)=closestBall(s(i).centroids,bLow(i-1,:,j));
% %     end
% % end

% % figure(1)
% % hold on
% % colorstring='brgmcykbrg';
% % for i=1:numberOfMarkersOnLowPlane
% %     plot(bLow(:,1,i),bLow(:,2,i),'Color',colorstring(i))
% % end
% % figure(2)
% % % imshow(im,'InitialMagnification','fit')
% % hold on
% % for i=1:numberOfMarkersOnLowPlane
% %     plot(bLow(:,1,i),bLow(:,2,i),'m')
% % end    

%*** original ***
% b1=zeros(N,2);
% b2=zeros(N,2);
% b1(1,:)=closestBall(s(1).centroids,[xtr(1),ytr(1)]);
% b2(1,:)=closestBall(s(1).centroids,[xtr(2),ytr(2)]);
% % b1(1,:)=s(1).centroids(4,:);
% % b2(1,:)=s(1).centroids(6,:);
% for i=2:N,
%     b1(i,:)=closestBall(s(i).centroids,b1(i-1,:));
%     b2(i,:)=closestBall(s(i).centroids,b2(i-1,:));
% end
% 
% figure(1)
% hold on
% plot(b1(:,1),b1(:,2),'r')
% plot(b2(:,1),b2(:,2),'c')
% 
% figure(2),clf
% imshow(im,'InitialMagnification','fit')
% hold on
% plot(b1(:,1),b1(:,2),'r')
% plot(b2(:,1),b2(:,2),'c')
%*** original ***

% return
%% Select 4 Reference points
%*** High Plane ***
% % % % % display('Select the 4 reference points for the High Plane,')
% % % % % figure(4),clf
% % % % % % im=imcrop(read(video,200),rect);
% % % % % im=imcrop(read(video,50),rect);
% % % % % imshow(im);
% % % % % % [xrefHigh,yrefHigh]=ginput(4);    % comment to use previous points (don't clear)
% % % % % % display([x,y])
% % % % % BW=detectColorHSV(im.*uint8(repmat(mask,[1,1,3])),4); % Yellow color
% % % % % clf
% % % % % imshow(BW)
% % % % % blob=regionprops(BW,'centroid');
% % % % % centroids=cat(1,blob.Centroid);
% % % % % ref1High=closestBall(centroids,[xrefHigh(1),yrefHigh(1)]);
% % % % % ref2High=closestBall(centroids,[xrefHigh(2),yrefHigh(2)]);
% % % % % ref3High=closestBall(centroids,[xrefHigh(3),yrefHigh(3)]);
% % % % % ref4High=closestBall(centroids,[xrefHigh(4),yrefHigh(4)]);
% % % % % figure(10),clf
% % % % % imshow(im)
% % % % % figure(11),clf
% % % % % imshow(BW)
%%
%*** Low Plane ***
% display('Select the 4 reference points for the Low Plane,')
% figure(4),clf
% % im=imcrop(read(video,200),rect);
% im=imcrop(read(video,50),rect);
% imshow(im);
% % [xrefLow,yrefLow]=ginput(4);    % comment to use previous points (don't clear)
% % display([x,y])
% BW=detectColorHSV(im.*uint8(repmat(mask,[1,1,3])),4); % Yellow color
% clf
% imshow(BW)
% blob=regionprops(BW,'centroid');
% centroids=cat(1,blob.Centroid);
% ref1Low=closestBall(centroids,[xrefLow(1),yrefLow(1)]);
% ref2Low=closestBall(centroids,[xrefLow(2),yrefLow(2)]);
% ref3Low=closestBall(centroids,[xrefLow(3),yrefLow(3)]);
% ref4Low=closestBall(centroids,[xrefLow(4),yrefLow(4)]);
% figure(12),clf
% imshow(im)
% figure(13),clf
% imshow(BW)
%% Move origin to lower left corner
%******** Dimensions of cropped image
dim=size(im1cropped);
xmax=dim(2);
ymax=dim(1);

%******** Flip coordinates to have the origin in the lower left corner
for i=1:numberOfMarkersOnHighPlane
    bbHigh(:,:,i)=bHigh(:,:,i);
end
% for i=1:numberOfMarkersOnLowPlane         % No Low Plane in 3LinkFish
%     bbLow(:,:,i)=bLow(:,:,i);
% end
bbHigh(:,2,:)=-bbHigh(:,2,:)+ymax;
% bbLow(:,2,:)=-bbLow(:,2,:)+ymax;
% % % ref1High(2)=-ref1High(2)+ymax;   % Already done in processing frames
% % % ref2High(2)=-ref2High(2)+ymax;
% % % ref3High(2)=-ref3High(2)+ymax;
% % % ref4High(2)=-ref4High(2)+ymax;
% ref1Low(2)=-ref1Low(2)+ymax;              % No Low Plane in 3LinkFish
% ref2Low(2)=-ref2Low(2)+ymax;
% ref3Low(2)=-ref3Low(2)+ymax;
% ref4Low(2)=-ref4Low(2)+ymax;

%*** original ***
% bb1=b1;
% bb1(:,2)=-b1(:,2)+ymax;
% bb2=b2;
% bb2(:,2)=-b2(:,2)+ymax;
% ref1High(2)=-ref1High(2)+ymax;
% ref2High(2)=-ref2High(2)+ymax;
% ref3High(2)=-ref3High(2)+ymax;
% ref4High(2)=-ref4High(2)+ymax;
%*** original ***

%% Projective Transformation (camera coordinates to inertial frame)
%Find Transformation
% T=createProjectiveTransformation([ref1(1);ref1(2);ref2(1);ref2(2);...
%     ref3(1);ref3(2);ref4(1);ref4(2)],...
%     2.54*[0;0;0;102;181.5;113.5;181.5;-4.5]);

%Find Transformation
%*** High Plane ***   (1. Bottom Left, 2. Top Left, 3. Top Right, 4. Bottom Right)
% % % % THigh=createProjectiveTransformation([ref1High(1);ref1High(2);ref2High(1);ref2High(2);...
% % % %     ref3High(1);ref3High(2);ref4High(1);ref4High(2)],...
% % % %     1*[6;0;6;117.95;232;118.75;232;0]);     %already done

%Find Transformation
%*** Low Plane ***
% TLow=createProjectiveTransformation([ref1Low(1);ref1Low(2);ref2Low(1);ref2Low(2);...
%     ref3Low(1);ref3Low(2);ref4Low(1);ref4Low(2)],...
%     1*[0;6;0;111.95;238;111.75;238;6]);

%Projection
%*** High Plane ***
bHighReal=zeros(N,2,numberOfMarkersOnHighPlane);
% bLowReal=zeros(N,2,numberOfMarkersOnLowPlane);
for i=1:numberOfMarkersOnHighPlane
    bHighReal(:,:,i)=projectCoordinates(bbHigh(:,:,i),THigh);
end
refReal=projectCoordinates([ref1High;ref2High;ref3High;ref4High;ref1High],THigh)
%*** Low Plane ***
% for i=1:numberOfMarkersOnLowPlane
%     bLowReal(:,:,i)=projectCoordinates(bbLow(:,:,i),TLow);
% end
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

%*** original ***
% breal1=projectCoordinates(bb1,T);
% breal2=projectCoordinates(bb2,T);
% 
% figure(5)
% plot(breal1(:,1),breal1(:,2))
% hold on
% plot(breal2(:,1),breal2(:,2),'r')
% plot((breal1(:,1)+breal2(:,1))/2,(breal1(:,2)+breal2(:,2))/2,'m')
% xlabel('x (cm)')
% ylabel('y (cm)')
% refReal=projectCoordinates([ref1High;ref2High;ref3High;ref4High;ref1High],T)
% plot(refReal(:,1),refReal(:,2))
%*** original ***

%% Plot fixed length between markers
diff=bHighReal(:,:,1)-bHighReal(:,:,2);
for i=1:length(diff), d(i)=norm(diff(i,:)); end
if(plots)
    figure(6),clf
    subplot(311)
    plot(d)
    hold on
    plot([0,length(diff)],[headDist,headDist],'k');
    plot([0,length(diff)],[mean(d),mean(d)],':k');
    legend('raw','real','mean')
    title('Head Link markers'' fixed length')
    xlabel('frame')
    ylabel('Length (cm)')
end

diff=bHighReal(:,:,3)-bHighReal(:,:,4);
for i=1:length(diff), d(i)=norm(diff(i,:)); end
if(plots)
%     figure(7),clf
    subplot(312)
    plot(d)
    hold on
    plot([0,length(diff)],[bodyDist,bodyDist],'k');
    plot([0,length(diff)],[mean(d),mean(d)],':k');
    legend('raw','real','mean')
    title('Body Link markers'' fixed length')
    xlabel('frame')
    ylabel('Length (cm)')
end

diff=bHighReal(:,:,5)-bHighReal(:,:,6);
for i=1:length(diff), d(i)=norm(diff(i,:)); end
if(plots)
%     figure(8),clf
    subplot(313)
    plot(d)
    hold on
    plot([0,length(diff)],[tailDist,tailDist],'k');
    plot([0,length(diff)],[mean(d),mean(d)],':k');
    legend('raw','real','mean')
    title('Tail Link markers'' fixed length')
    xlabel('frame')
    ylabel('Length (cm)')
end

% diff=breal1-breal2;
% for i=1:length(diff), d(i)=norm(diff(i,:)); end
% figure
% plot(d)
% hold on
% plot([0,length(diff)],[29.2,29.2],'k');
% plot([0,length(diff)],[mean(d),mean(d)],':k');
% legend('raw','real','mean')
% title('Rotor markers'' fixed length')
% xlabel('frame')
% ylabel('Length (cm)')

% save('SomeName','bHighReal','bLowReal','bHigh','bLow','bbHigh','bbLow','N','fps','im1cropped','imLastCropped','numberOfMarkersOnHighPlane','numberOfMarkersOnLowPlane','startFrame','THigh','TLow','ymax','xmax')

% return    
%% Save data...
%save('savedData\3LinkSwimmerSymmetric_Spring_f140_A60_S01.mat','fps','bHighReal','f','im1cropped','imLastCropped','bHigh')
save([folder,'savedData\',file1,file_middle,file2,'.mat'],'fps','bHighReal','f','im1cropped','imLastCropped','bHigh')
save(['savedData\',file1,file_middle,file2,'.mat'],'fps','bHighReal','f','im1cropped','imLastCropped','bHigh')
% PostProcess_3LinkEllipse    
figure(6)   % See errors...
set (figure(6), 'Units', 'normalized', 'Position', [0,0,1,1]);
% load handel
% sound(y,Fs)
figure(6), pause(1),
% saveAsPdf(6,['Processed_MarkersDistances\FixedMarkers_',file1,file_middle,file2])
saveAsPdf(6,[folder,'Processed_MarkersDistances\FixedMarkers_',file1,file_middle,file2])
pause(3)
