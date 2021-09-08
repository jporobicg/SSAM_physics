
function [numLayers, S] = get_numLayers(bgmfile, dlev);


dlev;
maxLayers = length(dlev);
[nbox,nface,bid,cent,b_area,verts,iface, botz, ibox] = read_boxes(bgmfile);

S = [];
botz;
for box = 1:nbox
    box
    botz(box)
    depth= 0;
    numLayers(box) = 0;
    if abs(botz(box)) > sum(dlev)
        botz(box) = sum(dlev) ;
    end
    for layer = 1:maxLayers
        if depth >= abs(botz(box))
            break;
        end
        numLayers(box) = layer;
        depth = depth + dlev(layer + 1)
    end
	S = sprintf('%s %d', S, numLayers(box));
end
S
botz;
