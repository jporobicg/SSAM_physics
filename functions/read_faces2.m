
function [nulr,nupt1,nupt2] = read_faces2(nbox, nface, bid,verts, iface, fname)

% Get the number of uniques faces  - this will be the array size.
nulr = repmat(-9,nface,2);
nupt1 = repmat(nan,nface,2);
nupt2 = repmat(nan,nface,2);

fid = fopen(fname,'r');


tin = fgetl(fid);
nf = 0;

while ischar(tin) 
    if ~isempty(findstr('# Face',tin))
        face = sscanf(tin,'# Face number %d');
        face = face + 1; % values start at 1 not 0 in matlab 
    elseif ~isempty(findstr('.p1',tin))
        nn = findstr('.p1',tin);
        nupt1(face,:) = sscanf(tin((nn+4):end),'%f%f');
    elseif ~isempty(findstr('.p2',tin))
        nn = findstr('.p2',tin);
        nupt2(face,:) = sscanf(tin((nn+4):end),'%f%f');
    elseif ~isempty(findstr('.lr',tin))
        nn = findstr('.lr',tin);
        nulr(face,:) = sscanf(tin((nn+4):end),'%d%d');
    end
    tin = fgetl(fid);
end;


fclose(fid);






