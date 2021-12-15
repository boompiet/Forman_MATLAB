function cmplx = musical_iso(cmplx,BO_1)
% -------------------------------------------------------------------------
% musical_iso.m
% -------------------------------------------------------------------------
% Purpose:      Compute musical isomorphism
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

nmax = 3*sum(cmplx(1).num(2).val);
n = 3*cmplx(1).num(1).val; m = cmplx(2).num(2).val;
iA = zeros(nmax,1); jA = zeros(nmax,1); sA = zeros(nmax,1); cnt = 0;
for i = 1:cmplx(1).num(1).val
    p = cmplx(1).num(2).val(i);
    iA(cnt+1:cnt+p)         = 3*(i-1)+1;
    iA(cnt+p+1:cnt+2*p)     = 3*(i-1)+2;
    iA(cnt+2*p+1:cnt+3*p)   = 3*(i-1)+3;
    
    jA(cnt+1:cnt+p)         = cmplx(1).bndop(2).indx(i,1:p);
    jA(cnt+p+1:cnt+2*p)     = cmplx(1).bndop(2).indx(i,1:p);
    jA(cnt+2*p+1:cnt+3*p)   = cmplx(1).bndop(2).indx(i,1:p);
    
    A = cmplx(2).vvol(cmplx(1).bndop(2).indx(i,1:p)).*...
        cmplx(2).dir(cmplx(1).bndop(2).indx(i,1:p),:); R = (A'*A)\A';
    
%     disp('---------')
%     disp(cmplx(1).bc(i,:))
%     disp(R)
    
    sA(cnt+1:cnt+p)         = R(1,:);
    sA(cnt+p+1:cnt+2*p)     = R(2,:);
    sA(cnt+2*p+1:cnt+3*p)   = R(3,:);
    
    cnt = cnt+3*p;
end
cmplx(2).sharp = sparse(iA,jA,sA,n,m,nmax);

cmplx(2).flat = kron(speye(cmplx(2).num(2).val),ones(1,3))*...
    spdiags(reshape(cmplx(2).vvol'.*cmplx(2).dir',3*cmplx(2).num(2).val,1),0,...
    3*cmplx(2).num(2).val,3*cmplx(2).num(2).val)*...
    kron(0.5.*abs(BO_1),speye(3));
