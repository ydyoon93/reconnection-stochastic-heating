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

for t=400
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
    
    figure(1)
    subaxis(1,3,1,'marginleft',.05,'marginbottom',.13)
    y = 201:401;
    x = 86:116;
    real_y = (y-Ymax/2-1)*.2;
    real_x = (x-Xmax/2-1)*.15;
    [real_x, real_y] = meshgrid(real_x,real_y);

    sx = [-2:0.15:-0.5 0.5:0.15:2 -2:0.15:-0.5 0.5:0.15:2];
    sy = [-20*ones(1,size(sx,2)/4) 20*ones(1,size(sx,2)/4)...
        20*ones(1,size(sx,2)/4) -20*ones(1,size(sx,2)/4)];
    h1 = contourf(real_x,real_y,phi(y,x),20,'linestyle','none');
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    hold on
    h2 = streamline(real_x,real_y,Qxe(y,x),Qye(y,x),sx,sy);
    set(h2,'color','r','linewidth',1.5)
    h3 = streamline(real_x,real_y,Bx(y,x),By(y,x),sx,sy);
    set(h3,'color','w','linewidth',1.5)
    plot(zeros(1,Ymax),linspace(-20,20,Ymax),'linestyle','--','color','m','linewidth',2)
    plot(linspace(-2,2,Xmax),zeros(1,Xmax),'linestyle','--','color','k','linewidth',2)
%     daspect([1 4 4])
    xlabel('$x/d_e$','interpreter','latex')
    ylabel('$y/d_e$','interpreter','latex')
    title('(a)','interpreter','latex')
    set(gca,'fontsize',15)
    
    subaxis(1,3,2,'marginright',0.05)
    y=301;
    Exdl = -cumsum(Ex(y,:))*.2;
    uzeBydl = -cumsum(uze(y,:).*By(y,:))*.2;
    plot(-20+20/200:40/200:20-20/200,Exdl-Exdl(101),'linewidth',2)
    hold on
    plot(-20:40/200:20-40/200,By(y,:).^2/2,'linewidth',2,'linestyle','--')
    plot(-20+20/200:40/200:20-20/200,uzeBydl-uzeBydl(101),'linewidth',2,'linestyle','-.')
    legend({'$-\int E_x dx+const.$','$B_y^2/2$','$-\int u_{ez}B_y dx+const.$'},'Location','southwest','interpreter','latex')
    plot(linspace(-12,-5,2),0.1*ones(1,2),'marker','>','markerindices',2,'linewidth',2,'color','g')
    plot(linspace(5,12,2),0.1*ones(1,2),'marker','<','markerindices',1,'linewidth',2,'color','g')
    xlabel('$x/d_e$','interpreter','latex')
    ylabel('Normalized $\phi$','interpreter','latex')
    xlim([-20,20])
    title('(b)','interpreter','latex')
    set(gca,'fontsize',15)
    
   
    subaxis(1,3,3,'marginright',.01)
    x=101;
    phi_eff = uze(:,x).^2/2+uye(:,x).^2;
    DelPhiEff=-(phi_eff(3:end)-phi_eff(1:end-2))/2/0.2;
    Eydl = -cumsum(Ey(:,x))*0.2;
    plot(-56+56/560:112/560:56-56/560,Eydl(21:end-20)+phi_eff(301)-Eydl(301),'linewidth',2)
    hold on
    plot(-56:112/560:56-112/560,phi_eff(21:end-20),'linewidth',2,'linestyle','--')
    legend({'$-\int E_y dy+const.$','$(u_{ey}^2+u_{e}^2)/2$'},'Location','southwest','interpreter','latex')
    plot(linspace(-30,-15,2),0.035*ones(1,2),'marker','<','markerindices',1,'linewidth',2,'color','g')
    plot(linspace(15,30,2),0.035*ones(1,2),'marker','>','markerindices',2,'linewidth',2,'color','g')
    xlabel('$y/d_e$','interpreter','latex')
    ylabel('Normalized $\phi$','interpreter','latex')
    xlim([-56,56])
    title('(c)','interpreter','latex')
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