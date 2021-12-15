% -------------------------------------------------------------------------
% setup_plotting.m
% -------------------------------------------------------------------------
% Purpose:      Setup indices for plotting purposes
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% sort arrays for plotting purposes
% -------------------------------------------------------------------------
disp('>> sort arrays for plotting purposes')

% Forman vertices
pindx = [1:vcmplx(1).num(1).val];
eindx = pindx(end) + [1:vcmplx(2).num(2).val];
findx = eindx(end) + [1:vcmplx(3).num(3).val];
vindx = findx(end) + [1:vcmplx(4).num(4).val];
[pz,pzi] = sort(z(pindx)); [ez,ezi] = sort(z(eindx));
[fz,fzi] = sort(z(findx)); [vz,vzi] = sort(z(vindx));
vindx = findx(end) + vzi; findx = eindx(end) + fzi;
eindx = pindx(end) + ezi; pindx = pzi;

% Forman edges
peindx = [1:sum(vcmplx(1).num(2).val)];
efindx = peindx(end) + [1:sum(vcmplx(2).num(3).val)];
fvindx = efindx(end) + [1:sum(vcmplx(3).num(4).val)];
z2 = fcmplx(2).bc(:,3); [pez,pezi] = sort(z2(peindx)); 
[efz,efzi] = sort(z2(efindx));[fvz,fvzi] = sort(z2(fvindx));
fvindx = efindx(end) + fvzi; efindx = peindx(end) + efzi; peindx = pezi;