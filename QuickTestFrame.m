
frame=285;
figure(50),clf,imshow(read(video,frame)),
figure(51),clf,imshow(detectColorHSV(read(video,frame),4))
figure(52),clf,imshow(bwareaopen(detectColorHSV(read(video,frame),4),P,CONN))

%%
figure(50),clf
analyzeHSV(read(video,frame))