% Author: Rodrigo Abrajan Guerrero
% 
% This function converts an RGB image into HSV space and returns a mask
% with the objects (pixels) that are of the color of interest. The color
% detection is done by defining a max and min threshold value for each of
% the HSV components of the color of interest. Right now it detects up to 4
% colors, Red, Green, Blue, Yellow.

function BW = detectColorHSV(Im,color)

nColors=length(color); % How many colors are we detecting, up to 4

for i=1:nColors         % set the threshold values for all req'd colors
    switch color(i)
        case 1
        % *** Red ***
%             hueThresholdLow=0.01;
%             hueThresholdHigh=.9;
%             saturationThresholdLow=0.5;
%             saturationThresholdHigh=0.9;
%             valueThresholdLow=0.4;
%             valueThresholdHigh=0.9;

%             hueThresholdLow=0.9;        %Using this before Rakshit's Pool
%             hueThresholdHigh=1;
%             saturationThresholdLow=0.575;    % Was 0.5, 0.575 good for most
%             saturationThresholdHigh=0.95;
%             valueThresholdLow=0.7;
%             valueThresholdHigh=1;
            
            %Trying for Rakshit's Pool (Red Ping-Pong Ball)
            hueThresholdLow=0.882;
            hueThresholdHigh=1;
            saturationThresholdLow=0.92;    % Was 0.5, 0.575 good for most
            saturationThresholdHigh=1;
            valueThresholdLow=0.17;
            valueThresholdHigh=0.45;
            
        case 2       
            % *** Green ***
            hueThresholdLow=0.1;
            hueThresholdHigh=.7;
            saturationThresholdLow=0.1;
            saturationThresholdHigh=0.9;
            valueThresholdLow=0;
            valueThresholdHigh=0.5;
        case 3
            % *** Blue ***
            hueThresholdLow=0.55;
            hueThresholdHigh=.7;
            saturationThresholdLow=0.6;
            saturationThresholdHigh=1;
            valueThresholdLow=0.7;
            valueThresholdHigh=1;
    %         hueThresholdLow=0.4;
    %         hueThresholdHigh=.85;
    %         saturationThresholdLow=0.45;
    %         saturationThresholdHigh=1;
    %         valueThresholdLow=0.55;
    %         valueThresholdHigh=1;
        case 4
            % *** Yellow ***
%             hueThresholdLow=0.08;
%             hueThresholdHigh=.145;
%             saturationThresholdLow=0.4;
%             saturationThresholdHigh=1;
%             valueThresholdLow=0.55;
%             valueThresholdHigh=1;        


%             hueThresholdLow=0.05; %0.05             % I was using this until 4/19/17
%             hueThresholdHigh=.15;  %0.15
%             saturationThresholdLow=0.35; %0.35
%             saturationThresholdHigh=1;
%             valueThresholdLow=0.55; %0.55
%             valueThresholdHigh=1;    
            
            hueThresholdLow=0.09; %0.05             % Changed to this on 4/19/17
            hueThresholdHigh=.15;  %0.15
            saturationThresholdLow=0.5; %0.35
            saturationThresholdHigh=1;
            valueThresholdLow=0.66; %0.55
            valueThresholdHigh=1;    
            
            %Trying for Rakshit's Pool (Yellow Marker) %Doesn't work great
%             hueThresholdLow=0.1020; %0.05
%             hueThresholdHigh=.2164;  %0.15
%             saturationThresholdLow=0.0224; %0.35
%             saturationThresholdHigh=0.4453;
%             valueThresholdLow=0.7175; %0.55
%             valueThresholdHigh=0.9075;                
        otherwise
            disp('Color not defined')
            return
    end
    lowBound(i,:)=[hueThresholdLow,saturationThresholdLow,valueThresholdLow];
    upperBound(i,:)=[hueThresholdHigh,saturationThresholdHigh,valueThresholdHigh];
end

HSV=rgb2hsv(Im);    % Convert to HSV color space
H=HSV(:,:,1);
S=HSV(:,:,2);
V=HSV(:,:,3);

if(nColors==1)
ind=1;
BW=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
% figure(10)
% subplot(221), imshow(mask1)
end

if(nColors==2)
ind=1;
mask1=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
ind=2;
mask2=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
BW=mask1|mask2;
% figure(10)
% subplot(221), imshow(mask1), subplot(222), imshow(mask2)
end

if(nColors==3)
ind=1;
mask1=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
ind=2;
mask2=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
ind=3;
mask3=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
BW=mask1|mask2|mask3;
% figure(10)
% subplot(221), imshow(mask1), subplot(222), imshow(mask2), subplot(223), imshow(mask3)
end

if(nColors==4)
ind=1;
mask1=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
ind=2;
mask2=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
ind=3;
mask3=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
ind=4;
mask4=((V >= lowBound(ind,3) & (V <= upperBound(ind,3)))&(S >= lowBound(ind,2) & (S <= upperBound(ind,2)))&(H >= lowBound(ind,1) & (H <= upperBound(ind,1))));
BW=mask1|mask2|mask3|mask4;
% figure(10)
% subplot(221), imshow(mask1), subplot(222), imshow(mask2), subplot(223), imshow(mask3), subplot(224), imshow(mask4)
end


% smallestAcceptableArea=1;
% BW = bwareaopen(BW, smallestAcceptableArea);
% structuringElement = strel('disk', 8);
% BW = imclose(BW, structuringElement);
% 
% BW=uint8(BW);
% maskedImageR = BW .* Im(:,:,1);
% maskedImageG = BW .* Im(:,:,2);
% maskedImageB = BW .* Im(:,:,3);
% maskedRGBImage = cat(3, maskedImageR, maskedImageG, maskedImageB);
% imshow(maskedRGBImage)




% *** Old Version, single color ***
% hueMask = (H >= hueThresholdLow) & (H <= hueThresholdHigh);
% saturationMask = (S >= saturationThresholdLow) & (S <= saturationThresholdHigh);
% valueMask = (V >= valueThresholdLow) & (V <= valueThresholdHigh);% BW=(hueMask & saturationMask & valueMask);
