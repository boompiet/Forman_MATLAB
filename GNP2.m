% -------------------------------------------------------------------------
% GNP2.m
% -------------------------------------------------------------------------
% Purpose:      
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% plot solution with higher coefficients in random cells
% -------------------------------------------------------------------------
close all; disp('>> plot solution with higher coefficients in random cells')

% extra array values needed
for i=1:4
    for j=1:4
        vcmplx(i).num(j).cval = cumsum(vcmplx(i).num(j).val);
        vcmplx(i).num(j).ival = vcmplx(i).num(j).cval - vcmplx(i).num(j).val + 1;
    end
end

range_percnt = [0:0.005:1]; num_MC = 1000;
data = zeros(length(range_percnt),num_MC,4);
for cnt = 1:num_MC
    cntp = 0; area = 0;
    k1=10^-10;k2=10^-10;k3=10^-10;
    k  = spdiags([k1*ones(sum(vcmplx(2).num(1).val),1);
        k2*ones(sum(vcmplx(3).num(2).val),1);
        k3*ones(sum(vcmplx(4).num(3).val),1)],0,fcmplx(2).num(2).val,fcmplx(2).num(2).val);
    indx = randperm(vcmplx(3).num(3).val);
    i1 = 1; i2 = 1; ig = fcmplx(1).bc(int_indx,3);
    for percnt = range_percnt
        cntp = cntp+1;
        disp([cntp,length(range_percnt),cnt,num_MC])
        
        % identify random faces for increased diffusion, as well as the
        % associated edge indices
        i1 = max(1,i2); i2 = round(percnt*vcmplx(3).num(3).val);
        GNPindxf = sort(indx(i1:i2));
        GNPindxe = unique(vcmplx(3).bndop(2).indx(GNPindxf,:));
        GNPindxe = sort(GNPindxe(GNPindxe>0));
        
        icnt = sum(vcmplx(2).num(1).val);
        for i=1:length(GNPindxf)
            area = area + vcmplx(3).vvol(GNPindxf(i));
            for j=vcmplx(3).num(2).ival(GNPindxf(i)):vcmplx(3).num(2).cval(GNPindxf(i))
                k(icnt+j,icnt+j) = 10^0; end; end
        
        for i=1:length(GNPindxe)
            for j=vcmplx(2).num(1).ival(GNPindxe(i)):vcmplx(2).num(1).cval(GNPindxe(i))
                k(j,j) = 10^0; end; end
        
        % setup system and solve problem
        A = fCBS_1 * k * fBO_1; b = - A*(z==1);
        t_old = zeros(fcmplx(1).num(1).val,1); t_old(z==1) = 1;
        t_new = A(int_indx,int_indx)\b(int_indx); t_old(int_indx) = t_new;
        
        flux = abs(k * fs2 * fBO_1 * t_old);

        zsf_v = zeros(size(fcmplx(1).bc(:,3)));  zsf_v(fcmplx(1).bc(:,3)<small)=1;
        zsf_e = abs(fBO_1)*zsf_v;  zsf_e(zsf_e>1)=0;
        zflux = zsf_e.*flux;
        
        osf_v = zeros(size(fcmplx(1).bc(:,3)));  osf_v(fcmplx(1).bc(:,3)>1-small)=1;
        osf_e = abs(fBO_1)*osf_v;  osf_e(osf_e>1)=0;
        oflux = osf_e.*flux;
        
        disp([sum(zflux),sum(oflux),abs(sum(zflux)-sum(oflux))])
        data(cntp,cnt,1:4) = [sum(zflux),sum(oflux),abs(sum(zflux)-sum(oflux)),area];
    end
end

data_mean = mean(data,2);  close all

figure; 
subplot(2,3,1);
semilogy(100.*range_percnt,data_mean(:,:,1),'*-b',100.*range_percnt,data_mean(:,:,2),'o--r')
xlabel('% of total faces in complex')
ylabel('Average effective conductivity')

subplot(2,3,2);
semilogy(data_mean(:,:,4),data_mean(:,:,1),'*-b',data_mean(:,:,4),data_mean(:,:,2),'o--r')
xlabel('Average cumulative area of GNP')
ylabel('Average effective conductivity')



