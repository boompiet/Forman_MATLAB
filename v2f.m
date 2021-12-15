function fcmplx = v2f(fcmplx,vcmplx)
% -------------------------------------------------------------------------
% v2f.m
% -------------------------------------------------------------------------
% Purpose:      Compute Forman's extended complex from Voronoi mesh
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Formann vertices
% -------------------------------------------------------------------------
fcmplx(1).num(1).val = vcmplx(1).num(1).val + vcmplx(2).num(2).val +...
    vcmplx(3).num(3).val + vcmplx(4).num(4).val;
fcmplx(1).bc = [vcmplx(1).bc; vcmplx(2).bc; vcmplx(3).bc; vcmplx(4).bc];
%fcmplx(1).cc = [vcmplx(1).cc; vcmplx(2).cc; vcmplx(3).cc; vcmplx(4).cc];

% -------------------------------------------------------------------------
% Formann edges - vertex connectivity
% -------------------------------------------------------------------------
fcmplx(2).num(2).val = sum(vcmplx(1).num(2).val) + ...
    sum(vcmplx(2).num(3).val) + sum(vcmplx(3).num(4).val);
fcmplx(2).num(1).val = 2*ones(fcmplx(2).num(2).val,1);
fcmplx(2).bndop(1).indx = zeros(fcmplx(2).num(2).val,2);
fcmplx(2).bndop(1).sgn = zeros(fcmplx(2).num(2).val,2);
cnt = 0; fcmplx(1).break(1) = vcmplx(1).num(1).val; 
for i=1:vcmplx(2).num(2).val
    for j=1:vcmplx(2).num(1).val(i)
        fcmplx(2).bndop(1).indx(cnt+j,:) = ...
            [vcmplx(2).bndop(1).indx(i,j), fcmplx(1).break(1) + i];
        fcmplx(2).bndop(1).sgn(cnt+j,:) = ...
            [vcmplx(2).bndop(1).sgn(i,j), -vcmplx(2).bndop(1).sgn(i,j)];
    end
    cnt = cnt + vcmplx(2).num(1).val(i);
end
fcmplx(2).break(1) = cnt;
fcmplx(1).break(2) = fcmplx(1).break(1) + vcmplx(2).num(2).val;
for i=1:vcmplx(3).num(3).val
    for j=1:vcmplx(3).num(2).val(i)
        fcmplx(2).bndop(1).indx(cnt+j,:) = [fcmplx(1).break(1) + ...
            vcmplx(3).bndop(2).indx(i,j), fcmplx(1).break(2) + i];
        fcmplx(2).bndop(1).sgn(cnt+j,:) = ...
            [vcmplx(3).bndop(2).sgn(i,j), -vcmplx(3).bndop(2).sgn(i,j)];
    end
    cnt = cnt + vcmplx(3).num(2).val(i);
end
fcmplx(2).break(2) = cnt;
fcmplx(1).break(3) = fcmplx(1).break(2) + vcmplx(3).num(3).val;
for i=1:vcmplx(4).num(4).val
    for j=1:vcmplx(4).num(3).val(i)
        fcmplx(2).bndop(1).indx(cnt+j,:) = [fcmplx(1).break(2) + ...
            vcmplx(4).bndop(3).indx(i,j), fcmplx(1).break(3) + i];
        fcmplx(2).bndop(1).sgn(cnt+j,:) = ...
            [vcmplx(4).bndop(3).sgn(i,j),-vcmplx(4).bndop(3).sgn(i,j)];
    end
    cnt = cnt + vcmplx(4).num(3).val(i);
end

% -------------------------------------------------------------------------
% Formann faces - edge connectivity
% -------------------------------------------------------------------------
fcmplx(3).num(3).val = sum(vcmplx(1).num(3).val) + ...
    sum(vcmplx(2).num(4).val);
fcmplx(3).num(2).val = 4*ones(fcmplx(3).num(3).val,1);
fcmplx(3).bndop(2).indx = zeros(fcmplx(3).num(3).val,4);
fcmplx(3).bndop(2).sgn = zeros(fcmplx(3).num(3).val,4);

cnt = 0;
ref0 = reshape(cumsum(reshape(abs(vcmplx(2).bndop(1).sgn)',...
    [1,numel(vcmplx(2).bndop(1).sgn)])),size(vcmplx(2).bndop(1).sgn'))';
ref1 = reshape(cumsum(reshape(abs(vcmplx(3).bndop(2).sgn)',...
    [1,numel(vcmplx(3).bndop(2).sgn)])),size(vcmplx(3).bndop(2).sgn'))';
ref2 = reshape(cumsum(reshape(abs(vcmplx(4).bndop(3).sgn)',...
    [1,numel(vcmplx(4).bndop(3).sgn)])),size(vcmplx(4).bndop(3).sgn'))';
for i=1:vcmplx(3).num(3).val
    for j=1:vcmplx(3).num(1).val(i)
        [r,c] = find(vcmplx(2).bndop(1).indx(vcmplx(3).bndop(2).indx...
            (i,1:vcmplx(3).num(2).val(i)),:)==...
            vcmplx(3).bndop(1).indx(i,j));
        
        fcmplx(3).bndop(2).indx(cnt+j,1:4) = [...
            ref0(vcmplx(3).bndop(2).indx(i,r(1)),c(1)),...
            ref0(vcmplx(3).bndop(2).indx(i,r(2)),c(2)),...
            fcmplx(2).break(1) + ref1(i,r(2)),...
            fcmplx(2).break(1) + ref1(i,r(1))];
        
        fcmplx(3).bndop(2).sgn(cnt+j,1:4) = [...
            vcmplx(3).bndop(2).sgn(i,r),...
            -fcmplx(2).bndop(1).sgn(fcmplx(3).bndop(2).indx(cnt+j,2),2),...
            fcmplx(2).bndop(1).sgn(fcmplx(3).bndop(2).indx(cnt+j,1),1)];
    end
    cnt = cnt + vcmplx(3).num(1).val(i);
end
fcmplx(3).break(1) = cnt;

for i=1:vcmplx(4).num(4).val
    for j=1:vcmplx(4).num(2).val(i)
        [r,c] = find(vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx...
            (i,1:vcmplx(4).num(3).val(i)),:)==...
            vcmplx(4).bndop(2).indx(i,j));
        
        f1 = fcmplx(2).break(1) + ...
            ref1(vcmplx(4).bndop(3).indx(i,r(1)),c(1));
        f2 = fcmplx(2).break(1) + ...
            ref1(vcmplx(4).bndop(3).indx(i,r(2)),c(2));
        f3 = fcmplx(2).break(2) + ref2(i,r(2));
        f4 = fcmplx(2).break(2) + ref2(i,r(1));
        
        fcmplx(3).bndop(2).indx(cnt+j,1:4) = [f1,f2,f3,f4];
        
        fcmplx(3).bndop(2).sgn(cnt+j,1:4) = [...
            fcmplx(2).bndop(1).sgn(f1,1),-fcmplx(2).bndop(1).sgn(f2,1),...
            -fcmplx(2).bndop(1).sgn(f3,1),fcmplx(2).bndop(1).sgn(f4,1)];
    end
    cnt = cnt + vcmplx(4).num(2).val(i);
end

% -------------------------------------------------------------------------
% Formann poly(hexa)hedra - face connectivity
% -------------------------------------------------------------------------
fcmplx(4).num(4).val = sum(vcmplx(4).num(1).val);
fcmplx(4).num(3).val = 6*ones(fcmplx(4).num(4).val,1);
fcmplx(4).bndop(3).indx = zeros(fcmplx(4).num(4).val,6);
fcmplx(4).bndop(3).sgn = zeros(fcmplx(4).num(4).val,6);

cnt = 0;
ref3 = reshape(cumsum(reshape(abs(vcmplx(3).bndop(1).sgn)',...
    [1,numel(vcmplx(3).bndop(1).sgn)])),size(vcmplx(3).bndop(1).sgn'))';
ref4 = reshape(cumsum(reshape(abs(vcmplx(4).bndop(2).sgn)',...
    [1,numel(vcmplx(4).bndop(2).sgn)])),size(vcmplx(4).bndop(2).sgn'))';
for i=1:vcmplx(4).num(4).val
    for j=1:vcmplx(4).num(1).val(i)
        [r,c] = find(vcmplx(3).bndop(1).indx(vcmplx(4).bndop(3).indx...
            (i,1:vcmplx(4).num(3).val(i)),:)==...
            vcmplx(4).bndop(1).indx(i,j));
        
        f1 = ref3(vcmplx(4).bndop(3).indx(i,r(1)),c(1));
        f2 = ref3(vcmplx(4).bndop(3).indx(i,r(2)),c(2));
        f3 = ref3(vcmplx(4).bndop(3).indx(i,r(3)),c(3));
        
        e4 = intersect(...
            vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx(i,r(1)),:),...
            vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx(i,r(2)),:));
        e5 = intersect(...
            vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx(i,r(2)),:),...
            vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx(i,r(3)),:));
        e6 = intersect(...
            vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx(i,r(3)),:),...
            vcmplx(3).bndop(2).indx(vcmplx(4).bndop(3).indx(i,r(1)),:));
        
        f4 = fcmplx(3).break(1) + ref4(i,...
            vcmplx(4).bndop(2).indx(i,1:vcmplx(4).num(2).val(i))==e4(end));
        f5 = fcmplx(3).break(1) + ref4(i,...
            vcmplx(4).bndop(2).indx(i,1:vcmplx(4).num(2).val(i))==e5(end));
        f6 = fcmplx(3).break(1) + ref4(i,...
            vcmplx(4).bndop(2).indx(i,1:vcmplx(4).num(2).val(i))==e6(end));
        
        fcmplx(4).bndop(3).indx(cnt+j,1:6) = [f1,f2,f3,f4,f5,f6];
        
        iref = fcmplx(3).bndop(2).indx([f1,f2,f3],3:4);
        sref = vcmplx(4).bndop(3).sgn(i,r)'.*...
            fcmplx(3).bndop(2).sgn([f1,f2,f3],3:4);
        
        fcmplx(4).bndop(3).sgn(cnt+j,1:6) = ...
            [vcmplx(4).bndop(3).sgn(i,r),...
            -fcmplx(3).bndop(2).sgn(f4,1)*...
            sref(iref==fcmplx(3).bndop(2).indx(f4,1)),...
            -fcmplx(3).bndop(2).sgn(f5,1)*...
            sref(iref==fcmplx(3).bndop(2).indx(f5,1)),...
            -fcmplx(3).bndop(2).sgn(f6,1)*...
            sref(iref==fcmplx(3).bndop(2).indx(f6,1))];
    end
    cnt = cnt + vcmplx(4).num(1).val(i);
end




