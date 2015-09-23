% generate artificial flow field rigid body flow
Rmax=30;
xvector=-Rmax:1:Rmax;
yvector=-Rmax:1:Rmax;
[xmesh,ymesh]=meshgrid(xvector,yvector);
omega=4;
data_u=-1*omega.*ymesh; %v=w*r
data_v=omega.*xmesh;
data_u=repmat(data_u,[1,1,100]);
data_v=repmat(data_v,[1,1,100]);

figure(1),
quiver(xmesh,ymesh,data_u(:,:,1),data_v(:,:,1),2);axis equal
title(' velocity field: rigid body,\omega=4 rad/s')
xlabel('x(mm)');
ylabel('y(mm)');

v_theta_temp=omega.*(1:Rmax);
v_theta=[-1*fliplr(v_theta_temp),0,v_theta_temp];
figure(2),
plot(xmesh(Rmax+1,:),v_theta,'-r');
hold on
plot(xmesh(Rmax+1,:),data_v(Rmax+1,:,1),'ob');
plot([-Rmax,Rmax],[0,0],'--k','linewidth',2.0)
hold off
title(' velocity comparison')
xlabel('x(mm)');
ylabel('tangential velocity(mm/s)');

%% caculate the acceleration
addpath('K:\PPIV')
delta_t=0.01;
[ax,ay,t_sel,xmesh_new,ymesh_new,xmesh_old,ymesh_old]=My_Get_Accelaration_centr_diff(xmesh,ymesh,data_u(:,:,1:3),data_v(:,:,1:3),delta_t,1,1,3);
figure(4),
quiver(xmesh_new,ymesh_new,ax,ay,2)


row_idx=Rmax+1;
aa=ax(row_idx,:,1);
bb=ay(row_idx,:,1);
aa2=-1*omega^2*xmesh(row_idx,:);
bb2=-1*omega^2*ymesh(row_idx,:);

figure(5),
plot(xmesh_new(row_idx,:),aa2,'-r','linewidth',2.0)
hold on
h(2)=plot(xmesh_new(row_idx,:),aa,'sk','markersize',8);
hold off
xlabel('r(m)');
ylabel('a_x(ms^{-2})');
legend('theory','caculated')

figure(6),
plot(xmesh_new(row_idx,:),bb2,'-r','linewidth',2.0)
hold on
h(2)=plot(xmesh_new(row_idx,:),bb,'ok','markersize',8);
hold off
xlabel('r(m)');
ylabel('a_y(ms^{-2})');
legend('theory','caculated')


%% caculate the pressure
tic
[pressure,xmesh_p,ymesh_p] = My_omni_directional_int(xmesh,ymesh,ax,ay,3);
toc

y=@(r)omega^2*r;%acceleration
pressure_theory_temp=zeros(Rmax,1);
for ii=1:Rmax
    pressure_theory_temp(ii)=quad(y,0,ii);
end
pressure_theory=[flipud(pressure_theory_temp);0;pressure_theory_temp];
 
figure(7),
pcolor(xmesh_p,ymesh_p,pressure);shading interp

figure(8),
pressure_caculated=pressure(row_idx-5,:)-min(min(pressure));
plot(-Rmax:Rmax,pressure_theory,'-r')
hold on
plot(xmesh_p(row_idx,:),pressure_caculated,'ob')
hold off
legend('theory','caculated')


