% Helper function to get the points in a repeating unit that are droplet

function points = get_droplet_points(unit_size, drop_radius)
    dropCentre = [(unit_size+1)/2, (unit_size+1)/2];
    dropletPoints = dropCentre;
    for i = 0 : floor(drop_radius)
        j = -floor(sqrt(drop_radius^2 - i^2)) : floor(sqrt(drop_radius^2 - i^2));
        dropletPoints = [dropletPoints; dropCentre + [repmat(i, [2*floor(sqrt(drop_radius^2 - i^2))+1, 1]), j']];
        if(i ~= 0)
            dropletPoints = [dropletPoints; dropCentre + [repmat(-i, [2*floor(sqrt(drop_radius^2 - i^2))+1, 1]), j']];
        end
    end

    points = unique(dropletPoints, 'rows');
end