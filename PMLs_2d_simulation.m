
nx = 100; 
ny = 100; 
nt = 500; 
dx = 10 / (nx - 1);
dy = 10 / (ny - 1);
dt = 0.01;
c = 3; 
u = zeros(nx, ny); 
v_x = zeros(nx+1, ny); 
v_y = zeros(nx, ny+1); 
psi = zeros(nx, ny); 


sigmax = zeros(nx+1, ny);
sigmay = zeros(nx, ny+1);

R = 1e-20;
m = 6    ;
source = zeros(nt,1);
sy  = floor(nx/2);
sx = floor(ny/2);
f0 = 2;   
t0 = 1.4 / f0;
for i = 1:nt
    tt = (i-1) * dt;
    source(i) = (1-2*(pi*f0*(tt-t0))^2)*exp(-(pi*f0*(tt-t0))^2);
end

ppw = c/(f0*dx);

for i = 1:nx
    if i >= 80
        value = -(c * log(R) * (m + 1)) / (20^(m + 1)) * (i - 80)^m;
        sigmax(i,:) = value;  
    elseif i <= 20
        value = -(c * log(R) * (m + 1)) / (20^(m + 1)) * (i - 20)^m;
        sigmax(i,:) = value;  
    end
end
for j = 1:ny
    if j >= 80
        value = -(c * log(R) * (m + 1)) / (20^(m + 1)) * (j - 80)^m;
        sigmay(:,j) = value;  
    elseif j <= 20
        value = -(c * log(R) * (m + 1)) / (20^(m + 1)) * (j - 20)^m;
        sigmay(:,j) = value;  
    end
end


filename = 'PMLs.gif';


for k = 1:nt
    for i = 2:nx
        for j = 1:ny
            v_x(i,j) = v_x(i,j) + dt * ((u(i,j) - u(i-1,j)) / dx - v_x(i,j) * sigmax(i,j));
        end
    end

    for i = 1:nx
        for j = 2:ny
            v_y(i,j) = v_y(i,j) + dt * ((u(i,j) - u(i,j-1)) / dy - v_y(i,j) * sigmay(i,j));
        end
    end

    for i = 2:nx-1
        for j = 2:ny-1
            psi(i,j) = psi(i,j) + dt * (sigmay(i,j) * (v_x(i+1,j) - v_x(i,j)) / dx + sigmax(i,j) * (v_y(i,j+1) - v_y(i,j)) / dy ...
                - (1 / c^2) * sigmax(i,j) * sigmay(i,j) * u(i,j));
        end
    end
    for i = 1:nx-1
        for j = 1:ny-1
            u(i,j) = u(i,j) + dt * (c^2 * ((v_x(i+1,j) - v_x(i,j)) / dx + (v_y(i,j+1) - v_y(i,j)) / dy ...
                - (1 / c^2) * (sigmax(i,j) + sigmay(i,j)) * u(i,j) + psi(i,j)));
            u(sy,sx) = u(sy,sx)+source(k);
        end
        
        
    end

    if mod(k, 10) == 0
        imagesc(u);  
        colorbar;
        title(['t = ', num2str(k * dt)]);
        xlabel('X');
        ylabel('Y');
        axis equal;
        clim([-1000,1000]);
        drawnow;
        
        % 捕捉图像帧并保存为GIF
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        
        if k == 10  % 第一帧
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1, 'DisposalMethod', 'restoreBG');
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1, 'DisposalMethod', 'restoreBG');
        end
    end
end
