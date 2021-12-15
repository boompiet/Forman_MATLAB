% -------------------------------------------------------------------------
% plot_sol_2D.m
% -------------------------------------------------------------------------
% Purpose:      Plot solution and flux in 2D plot
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
disp('>> plot solution 2D')

% solution 2D
figure;
subplot(2,2,1); plot(z(pindx),t_old(pindx),'k*',fcmplx(1).bc,fcmplx(1).bc,'r-');
xlabel('z-coordinate'); ylabel('Values'); title('Vertices'); legend('comp','exact','Location','northwest');
subplot(2,2,2); plot(z(eindx),t_old(eindx),'k*',fcmplx(1).bc,fcmplx(1).bc,'r-');
xlabel('z-coordinate'); ylabel('Values'); title('Edges'); legend('comp','exact','Location','northwest');
subplot(2,2,3); plot(z(findx),t_old(findx),'k*',fcmplx(1).bc,fcmplx(1).bc,'r-');
xlabel('z-coordinate'); ylabel('Values'); title('Faces'); legend('comp','exact','Location','northwest');
subplot(2,2,4); plot(z(vindx),t_old(vindx),'k*',fcmplx(1).bc,fcmplx(1).bc,'r-');
xlabel('z-coordinate'); ylabel('Values'); title('Volumes'); legend('comp','exact','Location','northwest');

% flux 2D
figure;
subplot(3,1,1); plot(z2(peindx),flux(peindx),'k*-');
xlabel('z-coordinate'); ylabel('Fluxes'); title('Vertices-Edges');
subplot(3,1,2); plot(z2(efindx),flux(efindx),'k*-');
xlabel('z-coordinate'); ylabel('Fluxes'); title('Edges-Faces');
subplot(3,1,3); plot(z2(fvindx),flux(fvindx),'k*-');
xlabel('z-coordinate'); ylabel('Fluxes'); title('Faces-Volumes');
