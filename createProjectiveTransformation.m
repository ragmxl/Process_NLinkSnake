function T=createProjectiveTransformation(x,xp)

% Author:Rodrigo Abrajan Guerrero
% 
% Used in conjuction with function projectCoordinates(posC,T), to give back the real
% coordinates from the coordinates obtained from a camera. All coordinates
% should be on the same plane.
% 
% 4 known positions are required to determine the transformation T.
% x represents the coordinates in the camera frame.
% xp are the coordinates of the 4 known positions in the inertial/real
% frame.


x1=x(1);
y1=x(2);
x2=x(3);
y2=x(4);
x3=x(5);
y3=x(6);
x4=x(7);
y4=x(8);

xp1=xp(1);
yp1=xp(2);
xp2=xp(3);
yp2=xp(4);
xp3=xp(5);
yp3=xp(6);
xp4=xp(7);
yp4=xp(8);

H=[ x1,y1,1,0,0,0,-x1*xp1,-y1*xp1;
    0,0,0,x1,y1,1,-x1*yp1,-y1*yp1;
    x2,y2,1,0,0,0,-x2*xp2,-y2*xp2;
    0,0,0,x2,y2,1,-x2*yp2,-y2*yp2;
    x3,y3,1,0,0,0,-x3*xp3,-y3*xp3;
    0,0,0,x3,y3,1,-x3*yp3,-y3*yp3;
    x4,y4,1,0,0,0,-x4*xp4,-y4*xp4;
    0,0,0,x4,y4,1,-x4*yp4,-y4*yp4 ];

a=inv(H)*xp;

T=[a(1),a(2),a(3);
   a(4),a(5),a(6);
   a(7),a(8),1];    