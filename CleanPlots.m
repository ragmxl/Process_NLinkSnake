clearvars -except ff
ff=8;
clc
file1='01.13.18\savedData\3LinkSnake_3Spring_A20_f';
file2='';
ext='.mat';
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
file=[file1,file_middle,file2,ext]
% return
load(file);
dt=1/fps;
t=0:dt:dt*(length(center)-1);
fontSize=18;

%%
figure(20),clf
plot(center(:,1),center(:,2))

% First move the initial position to the origin
points1=center;
points1=points1-[ones(length(points1),1)*center(1,1),ones(length(points1),1)*center(1,2)];
hold on
plot(points1(:,1),points1(:,2))

theta=atan2(b3(:,2)-b4(:,2),b3(:,1)-b4(:,1));
figure(21),clf
plot(theta*180/pi)

% Now, correct initial heading angle, so everything has theta = 0
R1=[cos(-theta(1)),-sin(-theta(1)); sin(-theta(1)),cos(-theta(1))];

xy=(R1*points1')';
figure(20)
plot(xy(:,1),xy(:,2)+10)
grid
axis equal

figure(22),clf
subplot(131)
plot(xy(:,1),xy(:,2)+10)
set(gca,'fontsize',fontSize)
grid
axis equal
xlabel('x (cm)')
ylabel('y (cm)')

subplot(132)
set(gca,'fontsize',fontSize)
plot(t,alpha1*180/pi)
set(gca,'fontsize',fontSize)
hold on
plot(t,alpha2*180/pi)
grid
ax=axis();
axis([0,1/f*5+1,ax(3),ax(4)])
legend('\alpha_1','\alpha_2')

subplot(133)
for cycle=1:5,
first=round(1+42/f*(cycle-1));
last=round(first+42/f);
% plot(alpha1,alpha2)
plot(alpha1(first:last)*180/pi,alpha2(first:last)*180/pi)
hold on
axis square
end
set(gca,'fontsize',fontSize)
grid on
legend('1T','2T','3T','4T','5T')
xlabel('\alpha_1')
ylabel('\alpha_2')
% xlabel('$\alpha_1$','Interpreter','LaTex')
% ylabel('$\alpha_2$','Interpreter','LaTex')






