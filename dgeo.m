function cmplx = dgeo(cmplx)
% -------------------------------------------------------------------------
% dgeo.m
% -------------------------------------------------------------------------
% Purpose:      Compute Delaunay geometry
%
% Pre:          Voronoi mesh has been computed
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% polyhedra - Delaunay vertices
% - cc (done previous)
% - dvol
% -------------------------------------------------------------------------
cmplx(4).dvol = ones(cmplx(4).num(4).val,1);

% -------------------------------------------------------------------------
% faces - Delaunay edges
% - cc
% - dvol
% -------------------------------------------------------------------------
cmplx(3).cc = zeros(cmplx(3).num(3).val,3);
cmplx(3).dvol = zeros(cmplx(3).num(3).val,1);

for i=1:cmplx(3).num(3).val
    % Delaunay edge circumcenter lies on Voronoi face.
    % Project Delaunay vertex along Voronoi face normal to face
    pnt = cmplx(4).cc(cmplx(3).bndop(4).indx(i,1),:);
    nrml = cmplx(3).bndop(4).sgn(i,1)*cmplx(3).dir(i,:);
    cmplx(3).cc(i,:) = pnt + dot(+cmplx(3).bc(i,:) - pnt,nrml)*nrml;

    % Delaunay edge length is the norm of the difference between end points
    pnt = cmplx(3).bndop(4).sgn(i,1)*pnt;
    if (cmplx(3).num(4).val(i)==2)
        pnt = pnt + cmplx(3).bndop(4).sgn(i,2)*...
            cmplx(4).cc(cmplx(3).bndop(4).indx(i,2),:);
    else; pnt = pnt - cmplx(3).bndop(4).sgn(i,1)*cmplx(3).cc(i,:);
    end
    cmplx(3).dvol(i) = norm(pnt);
end

% -------------------------------------------------------------------------
% edges - Delaunay faces
% - cc 
% - dvol
% -------------------------------------------------------------------------
cmplx(2).cc = zeros(cmplx(2).num(2).val,3);
cmplx(2).dvol = zeros(cmplx(2).num(2).val,1);

for i=1:cmplx(2).num(2).val
    % Projection of line from Delaunay edge circumcenter (on Delaunay face)
    % to Voronoi edge endpoint onto Voronoi edge
    epnt = cmplx(1).bc(cmplx(2).bndop(1).indx(i,1),:);
    ppnt = cmplx(3).cc(cmplx(2).bndop(3).indx(i,1),:);
    cmplx(2).cc(i,:) = epnt + ...
        dot((ppnt-epnt),cmplx(2).dir(i,:))*cmplx(2).dir(i,:);
    
    % area is the sum of triangles partitions (face,edge,vertex)
    % circumcenters
    for j=1:cmplx(2).num(3).val(i)
        k = cmplx(2).bndop(3).indx(i,j);
        d = norm(cmplx(2).cc(i,:) - cmplx(3).cc(k,:));
        sgn = sign(dot(cmplx(2).cc(i,:) - cmplx(3).cc(k,:),...
            cmplx(2).cc(i,:) - cmplx(3).bc(k,:)));
        cmplx(2).dvol(i) = cmplx(2).dvol(i) + sgn*0.5*d*cmplx(3).dvol(k);
        cmplx(2).dvol2(i,j) = sgn*0.5*d*cmplx(3).dvol(k);
    end
end

% -------------------------------------------------------------------------
% vertices - Delaunay tetra/hexahedra
% - cc (done previous; same as bc)
% - dvol
% -------------------------------------------------------------------------
cmplx(1).cc = cmplx(1).bc;
cmplx(1).dvol = zeros(cmplx(1).num(1).val,1);

for i=1:cmplx(1).num(1).val
    % volume is the sum of hexahedron partitions (tet,face,edge,vertex)
    % circumcenters
    for j=1:cmplx(1).num(2).val(i)
        k = cmplx(1).bndop(2).indx(i,j);
        d = norm(cmplx(1).cc(i,:) - cmplx(2).cc(k,:));
        sgn = sign(dot(cmplx(1).cc(i,:) - cmplx(2).cc(k,:),...
            cmplx(1).cc(i,:) - cmplx(2).bc(k,:)));
        cmplx(1).dvol(i) = cmplx(1).dvol(i) + sgn*d*cmplx(2).dvol(k)/3.0;
    end
end



