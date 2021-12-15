function BO = build_BO(cmplx,d)
% -------------------------------------------------------------------------
% build_BO.m
% -------------------------------------------------------------------------
% Purpose:      Build sparse boundary operator
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

n = sum(cmplx(d+1).num(d).val); m = cmplx(d+1).num(d+1).val;
iA = zeros(n,1); jA = zeros(n,1); sA = zeros(n,1); cnt = 0;
for i=1:m
    p = cmplx(d+1).num(d).val(i);
    iA(cnt+1:cnt+p) = i;
    jA(cnt+1:cnt+p) = cmplx(d+1).bndop(d).indx(i,1:p);
    sA(cnt+1:cnt+p) = cmplx(d+1).bndop(d).sgn(i,1:p);
    cnt = cnt+p;
end
BO= sparse(iA,jA,sA,m,cmplx(d).num(d).val,n);
