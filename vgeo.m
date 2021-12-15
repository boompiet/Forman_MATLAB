function cmplx = vgeo(cmplx)
% -------------------------------------------------------------------------
% vgeo.m
% -------------------------------------------------------------------------
% Purpose:      Compute Voronoi geometry
%
% Pre:          Neper mesh read using read_tess
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% nodes
% - bc (done previous)
% - dir (undefined)
% - vvol
% -------------------------------------------------------------------------
cmplx(1).vvol = ones(cmplx(1).num(1).val,1);

% -------------------------------------------------------------------------
% edges
% - bc 
% - dir
% - vvol
% -------------------------------------------------------------------------
cmplx(2).bc = 0.5*(cmplx(1).bc(cmplx(2).bndop(1).indx(:,1),:) + ...
    cmplx(1).bc(cmplx(2).bndop(1).indx(:,2),:));
cmplx(2).dir = ...
    cmplx(2).bndop(1).sgn(:,1).*...
    cmplx(1).bc(cmplx(2).bndop(1).indx(:,1),:) + ...
    cmplx(2).bndop(1).sgn(:,2).*...
    cmplx(1).bc(cmplx(2).bndop(1).indx(:,2),:);
cmplx(2).vvol = sqrt(sum((cmplx(2).dir).^2,2));
cmplx(2).dir = cmplx(2).dir./cmplx(2).vvol;

% -------------------------------------------------------------------------
% faces
% - bc
% - dir (done previous)
% - vvol
% -------------------------------------------------------------------------
cmplx(3).bc = zeros(cmplx(3).num(3).val,3);
cmplx(3).vvol = zeros(cmplx(3).num(3).val,1);
for i=1:cmplx(3).num(3).val
    for j=2:cmplx(3).num(1).val(i)-1
        tmpvol = 0.5*norm(cross(...
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,j),:) - ...
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,1),:),...
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,j+1),:)-...
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,1),:)));
        cmplx(3).bc(i,:) = cmplx(3).bc(i,:) + tmpvol*mean([...
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,1),:);
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,j),:);
            cmplx(1).bc(cmplx(3).bndop(1).indx(i,j+1),:)]);
        cmplx(3).vvol(i) = cmplx(3).vvol(i) + tmpvol;
    end
    cmplx(3).bc(i,:) = cmplx(3).bc(i,:)./cmplx(3).vvol(i);
end

% -------------------------------------------------------------------------
% polyhedra
% - bc
% - dir (undefined)
% - vvol
% -------------------------------------------------------------------------
cmplx(4).bc = zeros(cmplx(4).num(4).val,3);
cmplx(4).vvol = zeros(cmplx(4).num(4).val,1);
for i=1:cmplx(4).num(4).val
    refpnt = mean(cmplx(3).bc(cmplx(4).bndop(3).indx...
        (i,1:cmplx(4).num(3).val(i)),:),1);
    for j=1:cmplx(4).num(3).val(i)
        k = cmplx(4).bndop(3).indx(i,j);
        tmpvol = cmplx(4).bndop(3).sgn(i,j)/3.0*cmplx(3).vvol(k)*...
            dot(cmplx(3).dir(k,:),(cmplx(3).bc(k,:)-refpnt));
        cmplx(4).bc(i,:) = cmplx(4).bc(i,:) + tmpvol*(...
            cmplx(3).bc(k,:) + 0.25*(refpnt-cmplx(3).bc(k,:)));
        cmplx(4).vvol(i) = cmplx(4).vvol(i) + tmpvol;
    end
    cmplx(4).bc(i,:) = cmplx(4).bc(i,:)./cmplx(4).vvol(i);
end




