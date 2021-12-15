function cmplx = setup_complex(dim)
% -------------------------------------------------------------------------
% setup_complex.m
% -------------------------------------------------------------------------
% Purpose:      Setup data structure for complexes (probably no necessary)
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% name
% -------------------------------------------------------------------------
cmplx(dim).name='';

% -------------------------------------------------------------------------
% size
% -------------------------------------------------------------------------
cmplx(dim).num(dim).val=[]; % count of elements
cmplx(dim).break=[];        % 

% -------------------------------------------------------------------------
% geometry
% -------------------------------------------------------------------------
cmplx(dim).cc=[];       % circumcenter
cmplx(dim).bc=[];       % barycentric center

cmplx(dim).dir=[];      % direction

cmplx(dim).vvol=[];     % Voronoi volumne
cmplx(dim).dvol=[];     % Delaunay volume
cmplx(dim).bfvol=[];    % base Forman volume
cmplx(dim).fvol=[];     % Forman volume (inner product)

cmplx(dim).omega=[];    % curvature for boundary correction

% -------------------------------------------------------------------------
% connectivity
% -------------------------------------------------------------------------
cmplx(dim).bndop(dim).indx=[];  % boundary operator index
cmplx(dim).bndop(dim).sgn=[];   % boundary operator sign

cmplx(dim).flat=[];             % flat musical isomorphism
cmplx(dim).sharp=[];            % sharp musical isomorphism
