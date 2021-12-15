function cmplx = fgeo(cmplx)
% -------------------------------------------------------------------------
% fgeo.m
% -------------------------------------------------------------------------
% Purpose:      Compute Forman geometry
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% nodes
% - bfvol
% -------------------------------------------------------------------------
cmplx(1).bfvol = ones(cmplx(1).num(1).val,1);
cmplx(1).dfvol = ones(cmplx(1).num(1).val,1);

% -------------------------------------------------------------------------
% edges
% - dir
% - bfvol
% -------------------------------------------------------------------------
cmplx(2).bc = zeros(cmplx(2).num(2).val,3);
cmplx(2).dir = zeros(cmplx(2).num(2).val,3);
cmplx(2).bfvol = zeros(cmplx(2).num(2).val,1);
cmplx(2).dfvol = zeros(cmplx(2).num(2).val,1);
for i=1:cmplx(2).num(2).val
    cmplx(2).bc(i,:) = mean(cmplx(1).bc(cmplx(2).bndop(1).indx(i,1:2),:));
    
    cmplx(2).dir(i,:) = ...
        cmplx(2).bndop(1).sgn(i,1).*...
        cmplx(1).bc(cmplx(2).bndop(1).indx(i,1),:) + ...
        cmplx(2).bndop(1).sgn(i,2).*...
        cmplx(1).bc(cmplx(2).bndop(1).indx(i,2),:);
    
    cmplx(2).bfvol(i) = norm(cmplx(2).dir(i,:));
    cmplx(2).dfvol(i) = norm(cmplx(2).dir(i,:));
    
    cmplx(2).dir(i,:) = cmplx(2).dir(i,:) ./ cmplx(2).bfvol(i);
end

% -------------------------------------------------------------------------
% faces
% - dir
% - bfvol
% -------------------------------------------------------------------------
cmplx(3).bc = zeros(cmplx(3).num(3).val,3);
cmplx(3).dir = zeros(cmplx(3).num(3).val,3);
cmplx(3).bfvol = zeros(cmplx(3).num(3).val,1);
cmplx(3).dfvol = zeros(2*cmplx(3).num(3).val,1);
for i=1:cmplx(3).num(3).val
    cmplx(3).bc(i,:) = mean(cmplx(1).bc(cmplx(3).bndop(1).indx(i,1:4),:));
    
    for j=1:cmplx(3).num(2).val(i)
        k = cmplx(3).bndop(2).indx(i,j);
        
        h = norm(cross(cmplx(3).bc(i,:) - ...
            cmplx(2).bc(k,:),cmplx(2).dir(k,:)));
        cmplx(3).bfvol(i) = cmplx(3).bfvol(i) + 0.5*cmplx(2).bfvol(k)*h;
        
        if j==1 || j==4
            cmplx(3).dfvol(2*(i-1)+1) = ...
                cmplx(3).dfvol(2*(i-1)+1) + 0.5*cmplx(2).dfvol(k)*h;
        else
            cmplx(3).dfvol(2*(i-1)+2) = ...
                cmplx(3).dfvol(2*(i-1)+2) + 0.5*cmplx(2).dfvol(k)*h;
        end
    end
end

% -------------------------------------------------------------------------
% polyhedra
% - bfvol
% -------------------------------------------------------------------------
cmplx(4).bfvol = zeros(cmplx(4).num(4).val,1);
for i=1:cmplx(4).num(4).val
    cmplx(4).bc(i,:) = mean(cmplx(1).bc(cmplx(4).bndop(1).indx(i,:),:));
    for j=1:cmplx(4).num(3).val(i)
        k = cmplx(4).bndop(3).indx(i,j);
        for l=1:cmplx(3).num(2).val(k)
            m=cmplx(3).bndop(2).indx(k,l);
            edgs=[cmplx(3).bc(k,:);...
                cmplx(1).bc(cmplx(2).bndop(1).indx(m,1:2),:)] - ...
                cmplx(4).bc(i,:);
            cmplx(4).bfvol(i) = cmplx(4).bfvol(i) + ...
                cmplx(4).bndop(3).sgn(i,j)*cmplx(3).bndop(2).sgn(k,l)/6*...
                (dot(edgs(1,:),cmplx(2).bndop(1).sgn(m,1)*...
                cross(edgs(2,:),edgs(3,:))));
        end
    end
end
cmplx(4).bfvol=abs(cmplx(4).bfvol);

% -------------------------------------------------------------------------
% vertices
% - fvol
% -------------------------------------------------------------------------
cmplx(1).fvol  = zeros(cmplx(1).num(1).val,1);
cmplx(1).omega = zeros(cmplx(1).num(1).val,1);
for i=1:cmplx(1).num(1).val
    % boundary correction based on a curvature
    omega = 0;
    for j=1:cmplx(1).num(4).val(i)
        k = cmplx(1).bndop(4).indx(i,j); cnt = 1;
        v(1,1:3) = - cmplx(1).bc(i,:);
        v(2,1:3) = - cmplx(1).bc(i,:);
        v(3,1:3) = - cmplx(1).bc(i,:);
        for m=1:cmplx(4).num(2).val(k)
            n = cmplx(4).bndop(2).indx(k,m);
            if cmplx(2).bndop(1).indx(n,1)==i
                v(cnt,1:3) = v(cnt,1:3) + ...
                    cmplx(1).bc(cmplx(2).bndop(1).indx(n,2),:);
                cnt = cnt + 1;
            elseif cmplx(2).bndop(1).indx(n,2)==i
                v(cnt,1:3) = v(cnt,1:3) + ...
                    cmplx(1).bc(cmplx(2).bndop(1).indx(n,1),:);
                cnt = cnt + 1; end; end
        if cnt ~= 4; disp(cnt); pause; end
        a = v(1,1:3); b = v(2,1:3); c = v(3,1:3);
        tmp = 2*atan(abs(dot(a,cross(b,c))) / ...
            (norm(a)*norm(b)*norm(c) + ...
            dot(a,b)*norm(c) + dot(b,c)*norm(a) + dot(a,c)*norm(b)));
        if tmp<0; tmp = tmp + 2*pi; end
        omega = omega + tmp; cmplx(1).omega4(i,j) = tmp;
    end 
    % sum of volumes around nodes
    for j=1:cmplx(1).num(4).val(i)
        k = cmplx(1).bndop(4).indx(i,j);
        cmplx(1).fvol(i)  = cmplx(1).fvol(i)  + cmplx(4).bfvol(k); end
    % compute fvol
    omega = 4*pi()/omega;  cmplx(1).omega(i) = omega;
    cmplx(1).fvol(i)  = omega * cmplx(1).fvol(i)  / ...
        (2^(-0+3)*cmplx(1).bfvol(i)^2);
end

% -------------------------------------------------------------------------
% edges
% - fvol
% -------------------------------------------------------------------------
cmplx(2).fvol  = zeros(cmplx(2).num(2).val,1);
for i=1:cmplx(2).num(2).val
    
    for j=1:cmplx(2).num(1).val(i)
        k = cmplx(2).bndop(1).indx(i,j);
        
        cmplx(2).fvol(i) = cmplx(2).fvol(i) + cmplx(1).fvol(k);
    end
    cmplx(2).fvol(i)  = cmplx(2).fvol(i)  / (2^(1)*cmplx(2).bfvol(i)^2);
end

% -------------------------------------------------------------------------
% faces
% - fvol
% -------------------------------------------------------------------------
cmplx(3).fvol  = zeros(cmplx(3).num(3).val,1);
for i=1:cmplx(3).num(3).val
    for j=1:cmplx(3).num(1).val(i)
        k = cmplx(3).bndop(1).indx(i,j);
        
        cmplx(3).fvol(i) = cmplx(3).fvol(i) + cmplx(1).fvol(k);
    end
    cmplx(3).fvol(i)  = cmplx(3).fvol(i)  / (2^(2)*cmplx(3).bfvol(i)^2);
end

% -------------------------------------------------------------------------
% polyhedra
% - fvol
% -------------------------------------------------------------------------
cmplx(4).fvol  = zeros(cmplx(4).num(4).val,1);
for i=1:cmplx(4).num(4).val
    for j=1:cmplx(4).num(1).val(i)
        k = cmplx(4).bndop(1).indx(i,j);
        
        cmplx(4).fvol(i) = cmplx(4).fvol(i) + cmplx(1).fvol(k);
    end
    cmplx(4).fvol(i)  = cmplx(4).fvol(i)  / (2^(3)*cmplx(4).bfvol(i)^2);
end













