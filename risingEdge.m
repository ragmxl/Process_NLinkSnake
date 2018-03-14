function crossing = risingEdge(y,t)

% finds the times at which signal y crosses zero rising.
j=1;
for i=2:length(y),
    if(y(i)==0 && y(i-1)<0)     % Exactly ZERO
        crossing(j)=t(i);
        j=j+1;
    else if(y(i)>0 && y(i-1)<0)   % Zero crossing, find proportional time
            crossing(j)=(0-y(i))/(y(i)-y(i-1))*(t(i)-t(i-1))+t(i);
            j=j+1;            
        end
    end
end     