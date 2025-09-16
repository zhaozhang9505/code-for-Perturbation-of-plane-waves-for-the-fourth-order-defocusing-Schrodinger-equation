

% 孤子初值问题


clc
clear 
close all
x0=-100;
x1=50;
t=0:0.5:140;
N=128*2;

a=2;   %控制振幅

xx=linspace(x0,x1,N+1);
x=xx(1:N);
K=2*pi/(x1-x0)*[-N/2:N/2-1];
k=fftshift(K)';

n0=1-3/64*(4*cosh(1/8*x).^2+sinh(1/4*x).^2)./(2*cosh(1/4*x).*cosh(1/8*x)-sinh(1/4*x).*sinh(1/8*x)).^2;


v0=0.6778146454+3/64*(4*cosh(1/8*x).^2+sinh(1/4*x).^2)./(2*cosh(1/4*x).*cosh(1/8*x)-sinh(1/4*x).*sinh(1/8*x)).^2;



n0f=fft(n0);
v0f=fft(v0);

nvf=[n0f(:); v0f(:)];



[t,ufsol]=ode45('DNLSE4_zu',t,nvf,[],k,N);


nfsol=ufsol(:,1:N);

nsol=ifft(nfsol,[],2);


[X,T]=meshgrid(x,t);

subplot(1,2,1)
surf(X,T,-abs(nsol))
colormap(jet);
shading interp;
xlabel('\it \xi','FontSize',25,'FontWeight','bold');
ylabel('\it \tau','FontSize',25,'FontWeight','bold');
zlabel('-\it |q|^2','FontSize',25,'FontWeight','bold');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
view(41,41);
%colorbar
set(gcf,'color','w');
% zlim([-0.8 -1.05])
% clim([-0.9 -1])
 zlim([-1.05 -0.9])
 clim([-1 -0.9])
title('Numerical solution');





for i1=1:length(t)
    tau=t(i1);
    for i2=1:length(x)
        xi=x(i2);
        U1(i1,i2)=0.1e1 - 0.3e1 / 0.64e2 * (0.4e1 * cosh(xi / 0.8e1 + 0.1154601921e-1 * tau) ^ 2 + sinh(xi / 0.4e1 + 0.9236815366e-1 * tau) ^ 2) / (0.2e1 * cosh(xi / 0.4e1 + 0.9236815366e-1 * tau) * cosh(xi / 0.8e1 + 0.1154601921e-1 * tau) - sinh(xi / 0.4e1 + 0.9236815366e-1 * tau) * sinh(xi / 0.8e1 + 0.1154601921e-1 * tau)) ^ 2;
    end
end


subplot(1,2,2)
surf(X,T,-abs(U1));

colormap(jet);
clim([0 2]);
shading interp;
xlabel('\it \xi','FontSize',25,'FontWeight','bold');
ylabel('\it \tau','FontSize',25,'FontWeight','bold');
zlabel('-\it |q|^2','FontSize',25,'FontWeight','bold');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
view(41,41);
colorbar
set(gcf,'color','w');
% zlim([-0.8 -1.05])
% clim([-0.9 -1])
 zlim([-1.05 -0.9])
 clim([-1 -0.9])
set(gcf,'color','w');
title('Analytical solution');


%% 

figure

surf(X,T,abs(U1-nsol));

colormap(jet);
shading interp;
xlabel('\it x','FontSize',25,'FontWeight','bold');
ylabel('\it t','FontSize',25,'FontWeight','bold');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
colorbar
set(gcf,'color','w');
title('error');
xlim([-100 50])
zlim([0 0.01])
clim([0 0.01])


%% 
figure
z2=floor(length(t)/4);
z3=floor(length(t)/2);
z4=length(t);
z0=0.05;

subplot(2,2,1)
set(gcf,'color','w');
plot(X(1,:),abs(nsol(1,:)),'r',X(1,:),U1(1,:),'b','LineWidth',1.5)
title([' \tau = ' num2str(t(1))],'FontSize',15);
legend('Numerical solution','Analytical solution');
axis([min(x),max(x),min(min(abs(nsol)))-z0,max(max(abs(nsol)))+z0])
set(gca,'FontSize',15,'Fontname', 'Times New Roman')

subplot(2,2,2)
plot(X(z2,:),abs(nsol(z2,:)),'r',X(z2,:),U1(z2,:),'b','LineWidth',1.5)
title(['\tau = ' num2str(t(z2))],'FontSize',15);
legend('Numerical solution','Analytical solution');
axis([min(x),max(x),min(min(abs(nsol)))-z0,max(max(abs(nsol)))+z0])
set(gca,'FontSize',15,'Fontname', 'Times New Roman')

subplot(2,2,3)
plot(X(z3,:),abs(nsol(z3,:)),'r',X(z3,:),U1(z3,:),'b','LineWidth',1.5)
title(['\tau = ' num2str(t(z3))],'FontSize',15);
legend('Numerical solution','Analytical solution');
axis([min(x),max(x),min(min(abs(nsol)))-z0,max(max(abs(nsol)))+z0])
set(gca,'FontSize',15,'Fontname', 'Times New Roman')

subplot(2,2,4)
plot(X(z4,:),abs(nsol(z4,:)),'r',X(z4,:),U1(z4,:),'b','LineWidth',1.5)
title(['\tau = ' num2str(t(z4))],'FontSize',15);
legend('Numerical solution','Analytical solution');
axis([min(x),max(x),min(min(abs(nsol)))-z0,max(max(abs(nsol)))+z0])
set(gca,'FontSize',15,'Fontname', 'Times New Roman')









%% 
% plot(x,abs(nsol(1,:)))
% hold on
% plot(x,abs(nsol(1000,:)))
% plot(x,abs(nsol(5000,:)))


% p=polyfit(x,abs(nsol(8000,:)),120);
% figure
% zx=linspace(-150,40,800);
% zy=interp1(x,abs(nsol(24001,:)),zx,"cubic");
% %zy=polyval(p,zx);
% plot(zx,zy)


