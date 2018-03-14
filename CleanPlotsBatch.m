% This script plots the trajectories for a set of experiments, reorienting
% the trajectory such that the snake is alinged with the lab frame. For
% comparisson the starting positions are shifted vertically.

clc

figure(25),clf
figure(24),clf
figure(23),clf
figure(22),clf
figure(21),clf
figure(20),clf
leg={}
for ff=1:9,
    
    clearvars -except ff leg
%     ff=8;
    
%     file1='01.13.18\savedData\3LinkSnake_3Spring_A20_f';
%     title_='3 Link snake frequency sweep, A = 20 degrees, 3 "small" springs'; 
    
    file1='4Link\02.23.18\savedData\4LinkSnake_2_2_Spring_A20_f';
    title_='4 Link snake frequency sweep, A = 20 degrees, 2 & 2 springs'; 
    
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
    % First move the initial position to the origin
    points1=center;
    points1=points1-[ones(length(points1),1)*center(1,1),ones(length(points1),1)*center(1,2)];
%     hold on
%     plot(points1(:,1),points1(:,2))
    
    theta=atan2(b3(:,2)-b4(:,2),b3(:,1)-b4(:,1));
%     figure(21),clf
%     plot(theta*180/pi)
    
    % Now, correct initial heading angle, so everything has theta = 0
    R1=[cos(-theta(1)),-sin(-theta(1)); sin(-theta(1)),cos(-theta(1))];
    
    % Rotate 180 degrees to show start position in lower left corner
    R2=[cos(pi),-sin(pi); sin(pi),cos(pi)];
    
    xy=(R1*points1')';
    xy2=(R2*points1')';
 %%   
    figure(20),clf
    plot(xy2(:,1),xy2(:,2))
    set(gca,'fontsize',fontSize)
    title('Trejectory center of link 2')
    xlabel('x (cm)')
    ylabel('y (cm)')
    axis equal
    axis([0,230,0,130])
    grid
%     hold on
%     plot(xy2(1:5.15/f*fps,1),xy2(1:5.15/f*fps,2))
    figure(21),clf
    plot(t,alpha1*180/pi)
    set(gca,'fontsize',fontSize)
    hold on
    grid on
    plot(t,alpha2*180/pi)
    if(exist('alpha3'))
        plot(t,alpha3*180/pi)
        legend('\alpha_1','\alpha_2','\alpha_3')
    else
        legend('\alpha_1','\alpha_2')
    end
    xlabel('Time (s)')
    ylabel('Angle (degrees)')
    title('Joint angles')
    ax=axis;
    axis([ax(1),6/f,ax(3),ax(4)])
    figure(22),clf
    plot(alpha1*180/pi,alpha2*180/pi)
    set(gca,'fontsize',fontSize)
    xlabel('\alpha_1')
    ylabel('\alpha_2')
    title('Phase space joint angles')
    if(exist('alpha3'))
        figure(23),clf
        plot(alpha2*180/pi,alpha3*180/pi)
%         plot3(alpha1(3.5/f*fps:4.5/f*fps)*180/pi,alpha2(3.5/f*fps:4.5/f*fps)*180/pi,alpha3(3.5/f*fps:4.5/f*fps)*180/pi)
        set(gca,'fontsize',fontSize)
        xlabel('\alpha_2')
        ylabel('\alpha_3')
        title('Phase space joint angles')
                
        figure(24),clf
        plot3(alpha1*180/pi,alpha2*180/pi,alpha3*180/pi)
        set(gca,'fontsize',fontSize)
        xlabel('\alpha_1')
        ylabel('\alpha_2')
        zlabel('\alpha_3')
        title('Phase space joint angles')
        grid on
    else
        
    end
    
    figure(22)
%     keyboard
    
%%    
    
    figure(25)
    plot(xy(:,1),xy(:,2)+10*(ff-1))
    hold on
    set(gca,'fontsize',fontSize)
    grid
    axis equal
    xlabel('x (cm)')
    ylabel('y (cm)')
    leg(end+1)={[num2str(ff/10),' Hz']};    
end
grid on
legend(leg)
ax=axis;
axis([ax(1),ax(2),-10,10*ff])

% title('Frequency sweep, A=20 degrees, 3 "small" springs')
title(title_)





