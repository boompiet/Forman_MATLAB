% -------------------------------------------------------------------------
% plot_sol_2D.m
% -------------------------------------------------------------------------
% Purpose:      Plot solution and diffussion coefficient in 3D plot
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% plot solution
% -------------------------------------------------------------------------
disp('>> plot solution 3D')

% solution 3D
figure;
i = 1; j = vcmplx(1).num(1).val; subplot(2,2,1); set(gca,'FontSize',16);
plot3c2(fcmplx(1).bc(i:j,1),fcmplx(1).bc(i:j,2),fcmplx(1).bc(i:j,3),t_old(i:j),'o'); title('Vertices');
i = j+1; j = i-1+vcmplx(2).num(2).val; subplot(2,2,2); set(gca,'FontSize',16);
plot3c2(fcmplx(1).bc(i:j,1),fcmplx(1).bc(i:j,2),fcmplx(1).bc(i:j,3),t_old(i:j),'o'); title('Edges');
i = j+1; j = i-1+vcmplx(3).num(3).val; subplot(2,2,3); set(gca,'FontSize',16);
plot3c2(fcmplx(1).bc(i:j,1),fcmplx(1).bc(i:j,2),fcmplx(1).bc(i:j,3),t_old(i:j),'o'); title('Faces');
i = j+1; j = i-1+vcmplx(4).num(4).val; subplot(2,2,4); set(gca,'FontSize',16);
plot3c2(fcmplx(1).bc(i:j,1),fcmplx(1).bc(i:j,2),fcmplx(1).bc(i:j,3),t_old(i:j),'o'); title('Volumes');

% flux 3D
figure; kdiag = spdiags(k);
i = 1; j = sum(vcmplx(2).num(1).val); subplot(2,2,1); set(gca,'FontSize',16);
plot3c2(fcmplx(2).bc(i:j,1),fcmplx(2).bc(i:j,2),fcmplx(2).bc(i:j,3),kdiag(i:j),'o'); title('Vertices-Edges coefficients (k)');
i = j+1; j = i-1+sum(vcmplx(3).num(2).val); subplot(2,2,2); set(gca,'FontSize',16);
plot3c2(fcmplx(2).bc(i:j,1),fcmplx(2).bc(i:j,2),fcmplx(2).bc(i:j,3),kdiag(i:j),'o'); title('Edges-Faces coefficients (k)');
i = j+1; j = i-1+sum(vcmplx(4).num(3).val); subplot(2,2,3); set(gca,'FontSize',16);
plot3c2(fcmplx(2).bc(i:j,1),fcmplx(2).bc(i:j,2),fcmplx(2).bc(i:j,3),kdiag(i:j),'o'); title('Faces-Volumes coefficients (k)');

