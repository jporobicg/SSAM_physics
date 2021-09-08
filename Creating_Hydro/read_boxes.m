% read box definition bgm file
%
% jrd 13/10/05
%
% USAGE: [bid,cent,b_area,verts,iface] = read_boxes(fnm);

function [nbox, nface, bid,cent,b_area,verts,iface, b_botz, ibox] = read_boxes(fnm);

fid = fopen(fnm,'r');

verts{1} = [];
iface{1} = [];
b_area = [];
realarea = [];
b_botz = [];
cent = [];
bid = [];
ibox = [];
nbox = 0;
nface = 0;

nf = 0;

tin = fgetl(fid);

while ischar(tin) 
   fprintf('%s\n',tin);  

   if ~isempty(findstr('# Box number',tin))
      nf = nf+1;
      bid(nf) = sscanf(tin,'# Box number %d');
      nvert = 0;
   elseif ~isempty(findstr('nbox',tin))
      nbox = sscanf(tin,'nbox   %d');
   elseif ~isempty(findstr('nface',tin))
      nface = sscanf(tin,'nface %d');
   elseif ~isempty(findstr('.inside',tin))
      nn = findstr('.inside',tin);

      cent(nf,:) = sscanf(tin((nn+8):end),'%f%f');
   elseif ~isempty(findstr('.area',tin))
      nn = findstr('.area',tin);
      b_area(nf,:) = sscanf(tin((nn+6):end),'%f');
    elseif ~isempty(findstr('.botz',tin))
      nn = findstr('.botz',tin);
      b_botz(nf,:) = sscanf(tin((nn+6):end),'%f');
   elseif ~isempty(findstr('.vertmix',tin))
      % do nothing
   elseif ~isempty(findstr('.vert',tin))
      nvert = nvert+1;
      nn = findstr('.vert',tin);
      verts{nf}(nvert,:) = sscanf(tin((nn+6):end),'%f%f');
   elseif ~isempty(findstr('.iface',tin))
      nn = findstr('.iface',tin);
      iface{nf} = sscanf(tin((nn+7):end),'%d');
    elseif ~isempty(findstr('.ibox',tin))
      nn = findstr('.ibox',tin);
      ibox{nf} = sscanf(tin((nn+6):end),'%d') + 1;
   end
   tin = fgetl(fid);
end

%realarea = b_area.*cosd(cent(:,2));

fclose(fid);
