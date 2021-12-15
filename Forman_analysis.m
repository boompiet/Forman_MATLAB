% -------------------------------------------------------------------------
% Forman Analysis.m
% -------------------------------------------------------------------------
% Purpose:      Setup and solve diffusion equation using Forman's approach
%               as extended by Kiprian Berbatov
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% preliminaries
% -------------------------------------------------------------------------
clc;  close all;  clear all

% -------------------------------------------------------------------------
% create Voronoi complex data structure
% -------------------------------------------------------------------------
disp('>> create Voronoi complex structure')
vcmplx = setup_complex(4);
vcmplx(1).name='vertex'; vcmplx(2).name='edge';
vcmplx(3).name='face';   vcmplx(4).name='polyhedron';

% -------------------------------------------------------------------------
% read Neper .tess file for Voronoi complex
% -------------------------------------------------------------------------
filename = 'n100-id1.tess';                 % irregular mesh
% filename = 'nfrom_morpho-cube(2).tess';     % regular orthogonal mesh
disp(['>> read Neper .tess file (',filename,') for Voronoi complex'])
[vcmplx] = read_tess(filename,vcmplx);

% -------------------------------------------------------------------------
% get Voronoi geometry
% -------------------------------------------------------------------------
disp('>> get Voronoi geometry')
vcmplx = vgeo(vcmplx);

% -------------------------------------------------------------------------
% build Voronoi connectivity
% -------------------------------------------------------------------------
disp('>> build Voronoi connectivity')
vcmplx = vcon(vcmplx);
% vBO_1 = build_BO(vcmplx,1); % uncomment if you want the boundary
% vBO_2 = build_BO(vcmplx,2); % operators of the Voronoi complex
% vBO_3 = build_BO(vcmplx,3);

% -------------------------------------------------------------------------
% create Forman complex structure
% -------------------------------------------------------------------------
disp('>> create Forman complex structure')
fcmplx = setup_complex(4);
fcmplx(1).name='vertex'; fcmplx(2).name='edge';
fcmplx(3).name='face';   fcmplx(4).name='polyhedron';

% -------------------------------------------------------------------------
% setup Forman complex from Voronoi complex
% -------------------------------------------------------------------------
disp('>> setup Forman complex from Voronoi complex')
fcmplx = v2f(fcmplx,vcmplx);

% -------------------------------------------------------------------------
% build Forman connectivity - boundary operator fBO_1
% -------------------------------------------------------------------------
disp('>> build Forman connectivity')
fcmplx = fcon(fcmplx);
fBO_1 = build_BO(fcmplx,1);
% fBO_2 = build_BO(fcmplx,2); % uncomment if you need higher order
% fBO_3 = build_BO(fcmplx,3); % boundary operators

% -------------------------------------------------------------------------
% get Forman geometry
% -------------------------------------------------------------------------
disp('>> get Forman geometry')
fcmplx = fgeo(fcmplx);

% -------------------------------------------------------------------------
% get Forman coboundary star - fCBS_1
% -------------------------------------------------------------------------
disp('>> get Forman coboundary star')
innerNN_inv = spdiags(1./fcmplx(1).fvol,0,fcmplx(1).num(1).val,fcmplx(1).num(1).val);
innerEE  = spdiags(fcmplx(2).fvol,0,fcmplx(2).num(2).val,fcmplx(2).num(2).val);
fCBS_1 = - innerNN_inv * fBO_1' * innerEE;

% -------------------------------------------------------------------------
% set diffusion coefficients
% -------------------------------------------------------------------------
disp('>> set diffusion coefficients')
k1=1;k2=1;k3=1;
k  = spdiags([k1*ones(sum(vcmplx(2).num(1).val),1);
    k2*ones(sum(vcmplx(3).num(2).val),1);
    k3*ones(sum(vcmplx(4).num(3).val),1)],0,fcmplx(2).num(2).val,fcmplx(2).num(2).val);
% OPTIONAL: uncomment to use random diffusion coeffients
% k  = spdiags(10.^(5*([rand(sum(vcmplx(2).num(1).val),1);
%     rand(sum(vcmplx(3).num(2).val),1);
%     rand(sum(vcmplx(4).num(3).val),1)]-0.5)),0,fcmplx(2).num(2).val,fcmplx(2).num(2).val);

% -------------------------------------------------------------------------
% setup system and solve problem
%
% NOTE: this assumes a domain is a box of dimension 1x1x1
% -------------------------------------------------------------------------
disp('>> setup system and solve problem')
z  = fcmplx(1).bc(:,3); int_indx = find(z~=1 & z~=0);
A = fCBS_1 * k * fBO_1; b = - A*(z==1);
t_old = zeros(fcmplx(1).num(1).val,1); t_old(z==1) = 1;
t_new = A(int_indx,int_indx)\b(int_indx); t_old(int_indx) = t_new;

% -------------------------------------------------------------------------
% compute fluxes (effective diffusivity)
%
% NOTE: this assumes a domain is a box of dimension 1x1x1
% -------------------------------------------------------------------------
small=10^-8;
flux = abs(innerEE * k * fBO_1 * t_old);
zsf_v = zeros(size(fcmplx(1).bc(:,3)));  zsf_v(fcmplx(1).bc(:,3)<small)=1;  
zsf_e = abs(fBO_1)*zsf_v;  zsf_e(zsf_e>1)=0;  zflux = zsf_e.*flux;

osf_v = zeros(size(fcmplx(1).bc(:,3)));  osf_v(fcmplx(1).bc(:,3)>1-small)=1;  
osf_e = abs(fBO_1)*osf_v;  osf_e(osf_e>1)=0;  oflux = osf_e.*flux;

disp([sum(zflux),sum(oflux),abs(sum(zflux)-sum(oflux))])

% -------------------------------------------------------------------------
% plot solution
% -------------------------------------------------------------------------
setup_plotting;
plot_sol_2D;
plot_sol_3D;
return







