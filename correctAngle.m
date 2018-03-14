function y=correctAngle(phi)

y=zeros(length(phi),1);
y(1)=phi(1);
for i=2:length(phi),
    if(phi(i)-y(i-1) > pi/2)
        y(i)=phi(i)-pi;
    else if(phi(i)-y(i-1) < -pi/2)
            y(i)=phi(i)+pi;
        else
            y(i)=phi(i);
        end
    end
end