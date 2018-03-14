
% clear,clc
% file='savedData\3LinkSwimmerSymmetric_Spring_f150_A60_S01.mat';
% file='savedData\3LinkSwimmer_FatHead_Spring_f150_A60_S01.mat';
% file='savedData\3LinkSwimmer_SkinnyHead_Spring_f150_A60_S01.mat';
% file='savedData\3LinkSwimmer_FatHead_EqSpacing_Spring_f150_A60_S01.mat';

%% This block is used to find the initial frame, "startFrame"
% clearvars -except ff
clearvars -except amp

% ff=1;
clc
% 3LinkSnake_3Spring_A20_f01_012618.mat
file1='savedData/3LinkSnake_2Spring_A';
% file1='savedData\3LinkSwimmer_Symmetric_Spring_f';
% file2='_A60_S01';
% file2='_A60_S02';
file2='_f03';
ext='.mat';

ff=0.8;
f=ff/10             % Use ff in an external for loop

file_middle=num2str(amp);

% if(ff>=10)
%     file_middle=num2str(f*10);
% else if(ff>=1)
%         file_middle=['',num2str(f*10)];
%     else
%         file_middle=['0',num2str(f*100)];
%     end
% end

file=[file1,file_middle,file2,ext]
% file = 'savedData/3LinkSnake_3Spring_A20_f10.mat';
% return
load(file);

%%
startFrame=1;         % savedData includes this synced with experiment
pool_axes=[0,232,0,118.7];

figure(7),clf
imshow(im1cropped)
figure(8),clf
imshow(imLastCropped)

b1=bHighReal(startFrame:end,:,1);
b2=bHighReal(startFrame:end,:,2);
b3=bHighReal(startFrame:end,:,3);
b4=bHighReal(startFrame:end,:,4);
b5=bHighReal(startFrame:end,:,5);
b6=bHighReal(startFrame:end,:,6);

figure(9),clf
plot(b1(:,1),b1(:,2),'b')
hold on
plot(b2(:,1),b2(:,2),'r')
plot(b3(:,1),b3(:,2),'g')
plot(b4(:,1),b4(:,2),'m')
plot(b5(:,1),b5(:,2),'c')
plot(b6(:,1),b6(:,2),'y')
center=(b3+b4)/2;
plot(center(:,1),center(:,2),'k');
% WARNING: Axis values will need to change if lab frame coords change
axis([0 231.7 0 135.3]);
title('Marker Trajectories','Interpreter','latex','FontSize',20);
xlabel('x (cm)','Interpreter','latex','FontSize',16);
ylabel('y (cm)','Interpreter','latex','FontSize',16);
grid on
set(gca,'XTick',0:10:231.7);
set(gca,'YTick',0:10:135.3);

j=1;
for i=1:42:length(b1),       % every second...
    center_sample(j,:)=(b3(i,:)+b4(i,:))/2;
    j=j+1;
    i;
end
plot(center_sample(:,1),center_sample(:,2),'*k')
axis equal tight
axis(pool_axes)

% WARNING: Many of the values for resizing axes for display purposes will need to be
% changed if the laboratory reference frame changes
figure(20),clf
plot(center(:,1),center(:,2),'k');
axis([0 231.7 0 135.3]);
title('Centroid Trajectory','Interpreter','latex','FontSize',20);
xlabel('x (cm)','Interpreter','latex','FontSize',16);
ylabel('y (cm)','Interpreter','latex','FontSize',16);
grid on
set(gca,'XTick',0:10:231.7);
set(gca,'YTick',0:10:135.3);

%%
figure(10),clf
m_Head = (b1(:,2)-b2(:,2))./(b1(:,1)-b2(:,1));
% phi_Head2 = atan2(b1(:,2)-b2(:,2),b1(:,1)-b2(:,2));
phi_Head = atan(m_Head);
phi_Head = correctAngle(phi_Head);
m_Body = (b3(:,2)-b4(:,2))./(b3(:,1)-b4(:,1));
% phi_Body2 = atan2(b3(:,2)-b4(:,2),b3(:,1)-b4(:,1));
phi_Body = atan(m_Body);
m_Tail = (b5(:,2)-b6(:,2))./(b5(:,1)-b6(:,1));
phi_Tail = atan(m_Tail);
plot(phi_Head)
hold on
grid on
plot(phi_Body)
plot(phi_Head-phi_Body)
plot(phi_Tail)
plot(phi_Body-phi_Tail)
legend('Head','Body','HB','Tail','TB')
threshold=1*pi/180;
plot([1,400],[threshold,threshold],'--k')
plot([1,400],[-threshold,-threshold],'--k')

%%
figure(11),clf
dt=1/fps;
t=0:dt:dt*(length(b1)-1);
plot(t,phi_Head)
hold on
grid on
plot(t,phi_Body)
plot(t,phi_Head-phi_Body)
plot(t,phi_Tail)
plot(t,phi_Body-phi_Tail)
legend('Head','Body','HB','Tail','TB')
threshold=1*pi/180;
plot([0,t(end)],[threshold,threshold],'--k')
plot([0,t(end)],[-threshold,-threshold],'--k')

%% Plot Centroid Trajectory, stroboscopically sampled every once per cycle
%first find zero crossings in the joint angle...
jointAngle1=phi_Head-phi_Body;
times=risingEdge(jointAngle1,t);
strobed_center=interp1(t,center,times);
% j=0;
% T=1/f;
% for i=0:length(b1)-1,
%     if(t(i+1) >= j*T+0*T)  
%         strobed_center(j+1,:)=center(i+1,:);
%         j=j+1;
%     end
% end
figure(12),clf
plot(center(:,1),center(:,2))
hold on
plot(strobed_center(:,1),strobed_center(:,2),'*')
title(['Center position sampled with f = ',num2str(f),' Hz'],'Interpreter','latex','FontSize',20);
xlabel('x (cm)','Interpreter','latex','FontSize',16);
ylabel('y (cm)','Interpreter','latex','FontSize',16);
axis equal

%% Plot Joint Angle, stroboscopically sampled every once per cycle
strobed_jointAngle1=interp1(t,jointAngle1,times);
% j=0;
% T=1/f;
% for i=0:length(b1)-1,
%     if(t(i+1) >= j*T+0*T)  
%         strobed_jointAngle1(j+1,:)=jointAngle1(i+1,:);
%         strobed_t(j+1)=t(i+1);
%         j=j+1;
%     end
% end
figure(13),clf
plot(t,jointAngle1*180/pi)
hold on

plot(times,strobed_jointAngle1*180/pi,'*')
title(['Control input sampled with f = ',num2str(f),' Hz'],'Interpreter','latex','FontSize',20);
xlabel('t (s)','Interpreter','latex','FontSize',16);
ylabel('Amplitude (degrees)','Interpreter','latex','FontSize',16);
% axis equal
pause(1.5)

%% Plot joint angles in a separate figure
jointAngle2 = phi_Body-phi_Tail;

figure(14)
plot(t,jointAngle1*180/pi);
hold on
plot(t,jointAngle2*180/pi);
legend('Head Angle','Tail Angle');
title('Joint Angles','Interpreter','latex','FontSize',20);
xlabel('t (s)','Interpreter','latex','FontSize',16);
ylabel('Amplitude (degrees)','Interpreter','latex','FontSize',16);
axis([0 15 -amp-10 amp+10]);

figure(15),clf
%plot(jointAngle1*180/pi,jointAngle2*180/pi);
%plot(jointAngle1(100:end)*180/pi,jointAngle2(100:end)*180/pi)
plot(jointAngle1(250:540)*180/pi,jointAngle2(250:540)*180/pi)
title('Phase Space Plot','Interpreter','latex','FontSize',20)
xlabel('$$\alpha_1$$ (deg)','Interpreter','latex','FontSize',16);
ylabel('$$\alpha_2$$ (deg)','Interpreter','latex','FontSize',16);
pbaspect([1 1 1]);
grid on

%% To find the first frame "startFrame"...
return
figure(10)
axis([0,105,-0.2,0.2])
[x,y]=ginput(1);
startFrame=round(x)
% save(file,'-append','startFrame')     % Uncomment to save 'startFrame'

