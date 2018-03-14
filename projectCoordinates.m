function xReal=projectCoordinates(xCamera,T)

% Author: Rodrigo Abrajan Guerrero
% 
% Converts the coordinates given by a camera, xCamera, using the Projective
% Transformation T. 

% xpHomo=T*[xCamera(1); xCamera(2); 1];
% 
% xTemp=xpHomo/xpHomo(3);
% 
% xReal=[xTemp(1),xTemp(2)];
n=size(xCamera);
for i=1:n(1),
    xpHomo=T*[xCamera(i,1); xCamera(i,2); 1];
    xTemp=xpHomo/xpHomo(3);
%     xTemp=xpHomo;
    xReal(i,:)=[xTemp(1),xTemp(2)];
end
