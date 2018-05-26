Xmax = 200;
Ymax = 600;
Zmax = 50;
Tmax = 600;
mass_ratio = 10;
% 
% v = VideoWriter('C:\Users\ydyoo\Desktop\Reconnec_2D_py\Ion+Electron+Guide_Field_0+AR_1_1+MR_25\Bz_evolution.avi');
% open(v);
% scrsz = get(groot,'ScreenSize');
% figure('pos',[0 0 scrsz(3)/2 scrsz(4)*3/4])
% Qzi_fourier = zeros(Ymax,Tmax-1);
% Qze_fourier = zeros(Ymax,Tmax-1);
% Bz_fourier = zeros(Ymax,Tmax-1);
load('trajectory_2.mat')
for t=580
%     Qxi = data_read('Qxi.txt',Xmax,Ymax,t);
%     Qyi = data_read('Qyi.txt',Xmax,Ymax,t);
%     Qzi = data_read('Qzi.txt',Xmax,Ymax,t);
%     uxi = data_read('uxi.txt',Xmax,Ymax,t);
%     uyi = data_read('uyi.txt',Xmax,Ymax,t);
%     uzi = data_read('uzi.txt',Xmax,Ymax,t);
%     Qxe = data_read('Qxe.txt',Xmax,Ymax,t);
%     Qye = data_read('Qye.txt',Xmax,Ymax,t);
%     Qze = data_read('Qze.txt',Xmax,Ymax,t);
%     uxe = data_read('uxe.txt',Xmax,Ymax,t);
%     uye = data_read('uye.txt',Xmax,Ymax,t);
%     uze = data_read('uze.txt',Xmax,Ymax,t);
%     Bx = data_read('Bx.txt',Xmax,Ymax,t);
%     By = data_read('By.txt',Xmax,Ymax,t);
%     Bz = data_read('Bz.txt',Xmax,Ymax,t);
%     Ex = data_read('Ex.txt',Xmax,Ymax,t);
%     Ey = data_read('Ey.txt',Xmax,Ymax,t);
%     Ez = data_read('Ez.txt',Xmax,Ymax,t);
%     phi = data_read('phi.txt',Xmax,Ymax,t);
    
    y = 301:333;
    x = 89:113;
    real_y = (y-Ymax/2-1)*.2;
    real_x = (x-Xmax/2-1)*.15;
    figure(1);
    
    subaxis(2,4,1,'spacinghoriz',0.06,'spacingvert',0.12,'marginleft',0.03,'marginright',0.02)
    contourf(real_x,real_y,phi(y,x),100,'linestyle','none')
    hold on
    shift = 5700;
    plot(x_stable(:,(1:10)+shift),y_stable(:,(1:10)+shift),'linewidth',0.001,'color','r')
    plot(x_stable(:,1+shift),y_stable(:,1+shift),'linewidth',1.5,'color','w')
    xlabel('$x/d_e$','interpreter','latex')
    ylabel('$y/d_e$','interpreter','latex')
    title('(a)','interpreter','latex')
    set(gca,'fontsize',15)
    
    subaxis(2,4,5)
    contourf(real_x,real_y,phi(y,x),100,'linestyle','none')
    hold on
    plot(x_stochastic(:,(1:10)+shift),y_stochastic(:,(1:10)+shift),'linewidth',0.001,'color','r')
    plot(x_stochastic(:,1+shift),y_stochastic(:,1+shift),'linewidth',1.5,'color','w')
    xlabel('$x/d_e$','interpreter','latex')
    ylabel('$y/d_e$','interpreter','latex')
    title('(e)','interpreter','latex')
    set(gca,'fontsize',15)
    
    x1 = 3301;
    x2 = 3302;
    
    subaxis(2,4,2)
    scatter(y_stable(end,:),uy_stable(end,:),.7,'filled','markeredgecolor','k','markerfacecolor','k')
    hold on
    scatter(y_stable(end,x1),uy_stable(end,x1),30,'filled','markeredgecolor','r','markerfacecolor','r','markerfacealpha',1,'markeredgealpha',1)
    scatter(y_stable(end,x2),uy_stable(end,x2),30,'filled','markeredgecolor','r','markerfacecolor','r','markerfacealpha',1,'markeredgealpha',1)
    xlabel('$y/d_e$','interpreter','latex')
    ylabel('$\dot{y}/d_e|\omega|_{ce}$','interpreter','latex')
    xlim([.2 5])
    ylim([0 .05])
    title('(b)','interpreter','latex')
    set(gca,'fontsize',15)
    
    subaxis(2,4,6)
    scatter(y_stochastic(end,:),uy_stochastic(end,:),.7,'filled','markeredgecolor','k','markerfacecolor','k')
    hold on
    scatter(y_stochastic(end,x1),uy_stochastic(end,x1),30,'filled','markeredgecolor','r','markerfacecolor','r','markerfacealpha',1,'markeredgealpha',1)
    scatter(y_stochastic(end,x2),uy_stochastic(end,x2),30,'filled','markeredgecolor','r','markerfacecolor','r','markerfacealpha',1,'markeredgealpha',1)
    xlabel('$y/d_e$','interpreter','latex')
    ylabel('$\dot{y}/d_e|\omega|_{ce}$','interpreter','latex')
    xlim([.2 5])
    ylim([0 .05])
    title('(f)','interpreter','latex')
    set(gca,'fontsize',15)
    
    subaxis(2,4,3)
    scatter(y_stable(end,:),uy_stable(end,:),.7,'filled','markeredgecolor','k','markerfacecolor','k')
    xlabel('$y/d_e$','interpreter','latex')
    ylabel('$\dot{y}/d_e|\omega|_{ce}$','interpreter','latex')
    xlim([.2 1])
    ylim([0 .006])
    title('(c)','interpreter','latex')
    set(gca,'fontsize',15)
    
    subaxis(2,4,7)
    scatter(y_stochastic(end,:),uy_stochastic(end,:),.7,'filled','markeredgecolor','k','markerfacecolor','k')
    xlabel('$y/d_e$','interpreter','latex')
    ylabel('$\dot{y}/d_e|\omega|_{ce}$','interpreter','latex')
    xlim([.2 1])
    ylim([0 .006])
    title('(g)','interpreter','latex')
    set(gca,'fontsize',15)
    
    subaxis(2,4,8)
    histogram((uy_stable(end,:)-mean(uy_stable(end,:))).^2,20,'BinLimits',[0 .00004])
%     histogram((ux_stable(end,:)-mean(ux_stable(end,:))).^2,20)
    hold on
    histogram((uy_stochastic(end,:)-mean(uy_stochastic(end,:))).^2,20,'BinLimits',[0 .00004])
%     histogram((ux_stochastic(end,:)-mean(ux_stochastic(end,:))).^2,20)
    xlim([0 .00004])
    xlabel('$|\Delta\dot{y}|^2$','interpreter','latex')
    ylabel('Ion Counts','interpreter','latex')
    title('(h)','interpreter','latex')
    set(gca,'fontsize',15)
    
    
    subaxis(2,4,4)
    plot(sqrt((x_stable(:,x1)-x_stable(:,x2)).^2 ...
    +(y_stable(:,x1)-y_stable(:,x2)).^2 ...
    +(z_stable(:,x1)-z_stable(:,x2)).^2),'linewidth',3);
    hold on
    plot(sqrt((x_stochastic(:,x1)-x_stochastic(:,x2)).^2 ...
    +(y_stochastic(:,x1)-y_stochastic(:,x2)).^2 ...
    +(z_stochastic(:,x1)-z_stochastic(:,x2)).^2),'linewidth',3);
    xlabel('$t/|\omega_{ce}|^{-1}$','interpreter','latex')
    ylabel('Separation/$d_e$','interpreter','latex')
    title('(d)','interpreter','latex')
    set(gca,'fontsize',15)

%     Qxi_extend = replicate_in_3D(25*Qxi,Xmax,Ymax,Zmax);
%     Qyi_extend = replicate_in_3D(25*Qyi,Xmax,Ymax,Zmax);
%     Qzi_extend = replicate_in_3D(25*Qzi,Xmax,Ymax,Zmax);
% 
%     Qxe_extend = replicate_in_3D(Qxe,Xmax,Ymax,Zmax);
%     Qye_extend = replicate_in_3D(Qye,Xmax,Ymax,Zmax);
%     Qze_extend = replicate_in_3D(Qze,Xmax,Ymax,Zmax);
% 
%     Bx_extend = replicate_in_3D(Bx,Xmax,Ymax,Zmax);
%     By_extend = replicate_in_3D(By,Xmax,Ymax,Zmax);
%     Bz_extend = replicate_in_3D(Bz,Xmax,Ymax,Zmax);

%     stream_i = stream3(Qxi_extend,Qyi_extend,Qzi_extend,80,1,25);
%     stream_b = stream3(Bx_extend,By_extend,Bz_extend,90,1,25);
%     stream_e = stream3(Qxe_extend,Qye_extend,Qze_extend,80,600,25);
% 
%     si = streamtube(stream_i,5,[0.2 20]);
%     hold on
%     se = streamtube(stream_e,5,[0.2 20]);
%     sb = streamtube(stream_b,5,[0.2 20]);

    %% change the color of streamtubes

%     Qxi_init = data_read('Qxi.txt',Xmax,Ymax,1);
%     Qyi_init = data_read('Qyi.txt',Xmax,Ymax,1);
%     Qzi_init = data_read('Qzi.txt',Xmax,Ymax,1);
%     Qi_magnitude = sqrt(Qxi.^2+Qyi.^2+Qzi.^2)-sqrt(Qxi_init.^2+Qyi_init.^2+Qzi_init.^2);
%     Qi_magnitude = replicate_in_3D(Qi_magnitude,Xmax,Ymax,Zmax);
% 
%     for i=1:length(si)
%         %// Modify the colour data of each tube
%         set(si(i),'CData',interp3(Qi_magnitude,get(si(i),'XData')...
%             ,get(si(i),'YData'),get(si(i),'ZData'),'spline'))
%     end

%     Qxe_init = data_read('Qxe.txt',Xmax,Ymax,1);
%     Qye_init = data_read('Qye.txt',Xmax,Ymax,1);
%     Qze_init = data_read('Qze.txt',Xmax,Ymax,1);
%     Qe_magnitude = sqrt(Qxe.^2+Qye.^2+Qze.^2)-sqrt(Qxe_init.^2+Qye_init.^2+Qze_init.^2);
%     Qe_magnitude = replicate_in_3D(Qe_magnitude,Xmax,Ymax,Zmax);
% 
%     for i=1:length(se)
%         %// Modify the colour data of each tube
%         set(se(i),'CData',interp3(Qe_magnitude,get(se(i),'XData')...
%             ,get(se(i),'YData'),get(se(i),'ZData'),'spline'))
%     end


%     Bx_init = data_read('Bx.txt',Xmax,Ymax,1);
%     By_init = data_read('By.txt',Xmax,Ymax,1);
%     Bz_init = data_read('Bz.txt',Xmax,Ymax,1);
%     B_magnitude = sqrt(Bx.^2+By.^2+Bz.^2)-sqrt(Bx_init.^2+By_init.^2+Bz_init.^2);
%     B_magnitude = replicate_in_3D(B_magnitude,Xmax,Ymax,Zmax);
% 
%     for i=1:length(sb)
%         // Modify the colour data of each tube
%         set(sb(i),'CData',interp3(B_magnitude,get(sb(i),'XData')...
%             ,get(sb(i),'YData'),get(sb(i),'ZData'),'spline'))
%     end

%     caxis([-0.05 0.05])
%     colormap jet
%     colorbar
% 
%     shading interp
%     lighting gouraud
%     c1 = camlight;
%     xlim([77.5 81])
%     ylim([0 600])
%     zlim([17 31])
%     daspect([1 20 1])
%     grid on
%     view(3)

%     frame = getframe(gcf);
%     writeVideo(v,frame);

%     delete(si)
%     delete(se)
%     delete(sb)
%     delete(c1)
    disp(t)
end
% close(v)
% Qzi_fft = fftshift(abs(fft2(Qzi_fourier)));
% Qze_fft = fftshift(abs(fft2(Qze_fourier)));
% Bz_fft = fftshift(abs(fft2(Bz_fourier)));
%% curl u
% curl_uix = (Qxi-Bx)/mass_ratio;
% curl_uiy = (Qyi-By)/mass_ratio;
% curl_uiz = (Qzi-Bz)/mass_ratio;
% 
% curl_uex = Qxe+Bx;
% curl_uey = Qye+By;
% curl_uez = Qze+Bz;
% 
% curl_uix_extend = replicate_in_3D(curl_uix,Xmax,Ymax,Zmax);
% curl_uiy_extend = replicate_in_3D(curl_uiy,Xmax,Ymax,Zmax);
% curl_uiz_extend = replicate_in_3D(curl_uiz,Xmax,Ymax,Zmax);
% curl_uex_extend = replicate_in_3D(curl_uex,Xmax,Ymax,Zmax);
% curl_uey_extend = replicate_in_3D(curl_uey,Xmax,Ymax,Zmax);
% curl_uez_extend = replicate_in_3D(curl_uez,Xmax,Ymax,Zmax);
% 
% si = stream3(curl_uix_extend,curl_uiy_extend,curl_uiz_extend,83,1,30);
% sb = stream3(Bx_extend,By_extend,Bz_extend,83,1,30);
% se = stream3(curl_uex_extend,curl_uey_extend,curl_uez_extend,83,600,30);
% si = streamtube(si,5,[0.2 20]);
% hold on
% se = streamtube(se,5,[0.2 20]);
% sb = streamtube(sb,5,[0.2 20]);
% 
% shading interp
% lighting gouraud
% camlight
% xlim([77 82])
% ylim([0 600])
% zlim([15 32])
% daspect([1 20 1])
% grid on
% view(3)
%% for one-time run only
% figure(1)
% contourf(-0.5+2.5e-3:0.0025:0.5-2.5e-3,-2.5:(1/120):2.5-1/120,Qzi_fft,100,'linestyle','none')
% hold on
% x=0:0.0005:0.1;
% plot(x,sqrt(mass_ratio).*x,'color','r');
% plot(x,sqrt(x./(1-x)),'color','r');
% colorbar
% caxis([0 1])
% figure(2)
% contourf(-0.5+2.5e-3:0.0025:0.5-2.5e-3,-2.5:(1/120):2.5-1/120,Qze_fft,100,'linestyle','none')
% hold on
% plot(x,sqrt(mass_ratio).*x,'color','r');
% plot(x,sqrt(x./(1-x)),'color','r');
% colorbar
% caxis([0 1])
% figure(3)
% contourf(-0.5+2.5e-3:0.0025:0.5-2.5e-3,-2.5:(1/120):2.5-1/120,Bz_fft,100,'linestyle','none')
% hold on
% x=0:0.001:0.5;
% plot(x,sqrt(mass_ratio).*x,'color','r');
% plot(x,sqrt(x./(1-x)),'color','r');
% colorbar
% caxis([0 1])

%% functions
function data = data_read(filename,xmax,ymax,t)
    fileID = fopen(filename,'r');
    C = textscan(fileID,'%f',xmax*ymax,'Headerlines',t*xmax);
    data_temp = C{1,1};
    fclose(fileID);
    data = reshape(data_temp,[ymax xmax]);
end

function new_array = replicate_in_3D(array,xmax,ymax,zmax)
    new_array = zeros(ymax,xmax,zmax);
    for z = 1:zmax
        new_array(:,:,z) = array;
    end
end