function new = closestBall(list,old)


if(isempty(list)==0)
    new=list(1,:);
    closestDist=norm(old-new);
    dim=size(list);
    if(dim(1)>1)
        for i=2:length(list),
            tempDist=norm(old-list(i,:));
            if(tempDist<=closestDist)
                new=list(i,:);
                closestDist=tempDist;
            end        
        end
    end
else
    new=old;
end
    