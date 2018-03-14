% Author: Rodrigo Abrajan Guerrero
% 
% This function can be used to find the threshold values to do color
% detection using HSV space. This file was adapted from a script found
% online, "SimpleColorDetectionByHue.m" 


function analyzeHSV(im)

%select region to analyze
display('Select region to crop and analyze')
figure(40),clf
imshow(im)
[im,rect]=imcrop(im);

HSV=rgb2hsv(im);    % Convert to HSV color space
H=HSV(:,:,1);
S=HSV(:,:,2);
V=HSV(:,:,3);

figure(41),clf
h1=subplot(241);
imshow(im)
h2=subplot(242);
imshow(H)
h3=subplot(243);
imshow(S)
h4=subplot(244);
imshow(V)
linkaxes([h1,h2,h3,h4],'xy')

    fontSize=12;
    % Compute and plot the histogram of the "hue" band.
	hHuePlot = subplot(2, 4, 6); 
	[hueCounts, hueBinValues] = imhist(H); 
	maxHueBinValue = find(hueCounts > 0, 1, 'last'); 
	maxCountHue = max(hueCounts); 
	bar(hueBinValues, hueCounts, 'r'); 
	grid on; 
	xlabel('Hue Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Hue Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "saturation" band.
	hSaturationPlot = subplot(2, 4, 7); 
	[saturationCounts, saturationBinValues] = imhist(S); 
	maxSaturationBinValue = find(saturationCounts > 0, 1, 'last'); 
	maxCountSaturation = max(saturationCounts); 
	bar(saturationBinValues, saturationCounts, 'g', 'BarWidth', 0.95); 
	grid on; 
	xlabel('Saturation Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Saturation Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "value" band.
	hValuePlot = subplot(2, 4, 8); 
	[valueCounts, valueBinValues] = imhist(V); 
	maxValueBinValue = find(valueCounts > 0, 1, 'last'); 
	maxCountValue = max(valueCounts); 
	bar(valueBinValues, valueCounts, 'b'); 
	grid on; 
	xlabel('Value Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Value Image', 'FontSize', fontSize);
    
    	% Set all axes to be the same width and height.
	% This makes it easier to compare them.
	maxCount = max([maxCountHue,  maxCountSaturation, maxCountValue]); 
	axis([hHuePlot hSaturationPlot hValuePlot], [0 1 0 maxCount]); 

	% Plot all 3 histograms in one plot.
	subplot(2, 4, 5); 
	plot(hueBinValues, hueCounts, 'r', 'LineWidth', 2); 
	grid on; 
	xlabel('Values'); 
	ylabel('Pixel Count'); 
	hold on; 
	plot(saturationBinValues, saturationCounts, 'g', 'LineWidth', 2); 
	plot(valueBinValues, valueCounts, 'b', 'LineWidth', 2); 
	title('Histogram of All Bands', 'FontSize', fontSize); 
	maxGrayLevel = max([maxHueBinValue, maxSaturationBinValue, maxValueBinValue]); 
	% Make x-axis to just the max gray level on the bright end. 
	xlim([0 1]); 