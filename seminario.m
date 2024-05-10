
load('C:\Yop\Seminario\Threshold_2048_(current)\C2048_V1725dT30_coh.mat'); %reemplazar nombre en a

a =  'var_dataV1725dT30_incoh'; %cambiar el nombre al cambiar el archivo

a_ = 'varianzas_V1725dT30_incoh'; a_ = strrep(a_, '_', ' '); a_ = char(a_);


% calculo de varianzias en arreglos de 10x10  --> cada {i,j} 10x10 para todas las varianza de cada dot
centroidmat = TMs; 


M = 10;
var_datah = cell(M,M);
var_datav = cell(M,M);
%%
tic
for i=1:10
    for j=1:10
        var_datah{i, j} = wandering_struct(centroidmat,[i j],1);
        var_datav{i, j} = wandering_struct(centroidmat,[i j],2);
    end
end
toc
% guardando itvar, varianza total h + v para cada dot {i,j}
itvar = cell(M, M);
for i=1:10
    for j=1:10
        itvar{i,j} = var_datah{i,j} + var_datav{i,j};
    end
end
%
nombreh = [a '_horizontal'];
nombrev = [a '_vertical'];

eval([nombreh ' = var_datah;']);
eval([nombrev ' = var_datav;']);


%% referencias dots

centroid_mean = squeeze(mean(centroidmat,1));
ref_1 = reshape(centroid_mean',[],10, 10);
hor_refx = ref_1(1,1:10:end);hor_refy=ref_1(2,1:10:end);
ver_refx = ref_1(1,1:10); ver_refy = ref_1(2,1:10);
hor = polyfit(hor_refx,hor_refy,1);
ver = polyfit(ver_refx,ver_refy,1);

% angulos # especificar si es co o in 
angleh = atan(hor(1));anglev=atan(ver(1));
% Referencia  x e y  
REFX = squeeze(ref_1(1,:,:));REFY =squeeze(ref_1(2,:,:));
%% Calculo de aoa, BW 
aoa = aoa_struct(centroidmat,1)+ aoa_struct(centroidmat,2);
namevar = a;
rx = 'REFX'; ry = 'REFY'; wanvar = 'itvar'; aoavar = 'aoa'; 
eval(['rwvar' namevar '= struct(rx, REFX, ry, REFY, wanvar, itvar, aoavar, aoa);']);
%% vemos los puntos de referencia obtenidos de la medida
p = plot(REFX, REFY, '-s', 'MarkerSize', 10, 'MarkerFaceColor',[1 .6 .6], 'MarkerEdgeColor','blue');


%% mapeo de las variancias de posiciones

itvar40 = itvar{5,5};
figure; contourf(REFX,REFY,itvar40,'ShowText','on')
title(a_)

% hasta aqui tenemos las variancias de puntos de 1 medida



%% Grabar grafico - no correr, solo para ver animacion 
v = VideoWriter('C:\Yop\Seminario\animacion3D_varianzas_final.mp4', 'MPEG-4'); % Puedes elegir otros formatos y configuraciones según tus preferencias
% Configurar las propiedades del VideoWriter
v.FrameRate = 100; % Número de cuadros por segundo
open(v);
% grafico 3d
figure;
surf(REFX, REFY, itvar40); % representacion en 3d 
title('Representación 3D');
xlabel('Eje X');
ylabel('Eje Y');
zlabel('Magnitud de las varianzas');
colormap('hot')
view(45, 45);

% Animación en 2D (vista desde arriba)
for angulo = 0:1:360 % Cambia el ángulo según tu animación
    % Rotar la vista de los ejes en los ejes X e Y (vista desde arriba)
    view(0, 90); % Vista desde arriba
    
    % Capturar el cuadro actual de los ejes
    frame = getframe(gca);
    frame = frame.cdata;    
    frame = imresize(frame, [344, 434]);
    
    % Escribe el cuadro en el archivo de video
    writeVideo(v, frame);
end
% Realizar la animación
for angulo = 0:1:720 % Cambia el ángulo según tu animación
    % Rotar la vista de la figura 3D
    camorbit(0.1767, 0.1767, 'camera');
    
    % Capturar el cuadro actual y escribirlo en el video
    frame = getframe(gca);
    frame = frame.cdata;
    frame = imresize(frame, [344, 434]);

    writeVideo(v, frame);
end

% Cerrar el VideoWriter
close(v);



%% Parametros del experimento
% On the DMD
px_size = 20; 
size_E = 471; %max wide in px, each spot is 12x12 and the interspace is 39 px. between centers is
d_spot = 51; % is the separation between centers in the dots
N      = 10;  % number of dots

% On the image plane, and observed by the camera
MA     = 3.884; % is the magnification of the inverted telescope
esize_E= size_E*px_size*MA; % 2.5 cm at the plane
d_eff  = d_spot*px_size*MA; % 2.71mm at the plane


% Specifics of the image system
pL = 84;%cm is the distance between the camera and the image plane of the inverted telescope
f = 85;%mm is the focal for the camera
D = f/1.8; %47.22 mm according to the specifications in the lens.
pixel_scale  = 20e-6/.085;%   in radians. Angle of arrival is pixel_size/f.
sampling_freq = 1000;%Hz

% Supected parameters from the turbulence are
l0          = 0.008; %m is the inner scale (the fit should confirm or not these values)
lambda_incoh= 660e-9; %m for incoherent light targets
lambda_coh  = 637e-9; %m for coherent light
L           = 0.37;   %m is the turbulent chamber length
% To use the functions defined to estimate the variances between the LEDs for the incoherence case
Fm_incoh = pi*L*lambda_incoh/l0^2; %0.01102 []
Fm_coh   = pi*L*lambda_coh/l0^2;   %0.01060 []

Rm = pi/l0*D; %18.54
delta_m = 2*pi/l0*d_eff; %2.12, delta min.
Delta_m = 2*pi/l0*esize_E;%19.63 delta max

%% comparando

Deltax0= REFX(:,end)-REFX(:,1);Deltay0 = REFY(:,end)-REFY(:,1);
Deltax1= REFX(end,:)-REFX(1,:);Deltay1 = REFY(end,:)-REFY(1,:);
size_targeth = mean(sqrt(Deltax0.^2+Deltay0.^2),"all");
size_targetv = mean(sqrt(Deltax1.^2+Deltay1.^2),"all");
deltax0 = REFX(:,2:end)-REFX(:,1:end-1);deltay0 = REFY(:,2:end)-REFY(:,1:end-1);
deltax1 = REFX(2:end,:)-REFX(1:end-1,:);deltay1 = REFY(2:end,:)-REFY(1:end-1,:);
delta_targeth_co = mean(sqrt(deltax0.^2+deltay0.^2),"all")/2+size_targeth/18;
delta_targetv_co = mean(sqrt(deltax1.^2+deltay1.^2),"all")/2+size_targetv/18;
% solo por comparar y0=REFY(1,1);x0=REFX(1,1); 
clear Delta* deltax* deltay1

%% La estimación de lineas ya la tenemos
%% Angulos y delta_target 
% coorrection for real target sizes, because of obliquity

% angleh y anglev ya los obtuvimos, ...

%incoherente
%theta_in =anglev-angleh;
%delta_targetm2 = (delta_targetv_in^2+delta_targeth_in^2)/2;
%a_in = 9*sqrt(delta_targetm2+sqrt(delta_targetm2^2-delta_targetv_in^2*delta_targeth_in^2*sin(theta_in)^2));


%coherente - nuestra medida es coherente
theta_co =anglev-angleh;
delta_targetm2 = (delta_targetv_co^2+delta_targeth_co^2)/2;
a_co = 9*sqrt(delta_targetm2+sqrt(delta_targetm2^2-delta_targetv_co^2*delta_targeth_co^2*sin(theta_co)^2));

clear delta_targetv_m2 angleh anglev

% Correct zoom magnification
%ma = (px_size*MA/a_in+px_size*MA/a_co)*9*51/40;
ma = .1/(504.1676*20*10^(-6)); %magnificación
% Parameters for the lab.

%% Valores para calcular DA4 
% On the object plane parameters
px_dmd = 13.67; % µm
px_size = 20; % µm, at the CMOS
size_E = 471; %max wide in px, each spot is 12x12 and the interspace is 39 px. between centers is
d_spot = 51; % is the separation between centers in the dots
M      = 10;  % number of dots
% Supected parameters from the turbulence are
l0          = 0.008; %m is the inner scale (the fit should confirm or not these values)
L0          = 0.15;  %m is the outer scale
H           = 1/3;   % is the Kolmogorov exponent. We should try with the one estimated from mDFA
lambda_in= 660e-9; %m for incoherent light targets
lambda_co  = 637e-9; %m for coherent light
L           = 84;   %m is the turbulent chamber length
% On the image plane, and observed by the camera
Fm_co = pi*lambda_co/l0^2;Fm_in = pi*lambda_in/l0^2; %Fresnel Number for either incoherent or coherent
%deltamh_in=2*pi*delta_targeth_in*ma*px_size*10^(-6)/l0; %deltam horizontal
%deltamv_in=2*pi*delta_targetv_in*ma*px_size*10^(-6)/l0; %deltam vertical
deltamh_co=2*pi*delta_targeth_co*ma*px_size*10^(-6)/l0; %deltam horizontal
deltamv_co=2*pi*delta_targetv_co*ma*px_size*10^(-6)/l0; %deltam vertical
qlab = l0/L0;
D = 0.04722;
Rmlab = pi*D/l0;%47.22
kmlab = 2*pi/l0;
%pixel_scale  = 20e-6/.085;%   in radians. Angle of arrival is pixel_size/f.


%% Se crean sets de fits
% [est,fixrefx,fixrefy] = Sinv(nvar, range,row,col,hv, M)
mlist = ['08';'12';'17';'25';'30';'35';'40'];

% Set data to estimate for a given measure.

% nvar=['*V' mlist(7,:) mlist(7,:) 'dT*_co*']; %cambiar segun mi workspace
rwvar = 'rwvar';
nvar = [rwvar a];
xydatah_co = cell(2,M*M);
xydatav_co = cell(2,M*M);
for k=1:100
    i = floor((k - 1) / M) + 1;
    j = mod(k - 1, M) + 1;
    [est,~,~] = Sinv(nvar, 5,i,j,1, 10);
    [~,col] = find(est==0); % finds zeros
    est1 = fliplr(est(:,1:col(1))); % set straight first part
    est2 = est(:,col(1):end);
    if k==6
    disp(k);
    end
    sz1 = size(est1)-1;sz2 = size(est2)-1;
    ref_est1 = (0:sz1(2))*delta_targeth_co;
    ref_est2 = (0:sz2(2))*delta_targeth_co;
    est_co = [est1(:)' est2(:)'];
    ref_est1 = repelem(ref_est1,sz1(1)+1,1);
    ref_est2 = repelem(ref_est2,sz1(1)+1,1);
    ref_co = [ref_est1(:)',ref_est2(:)'];
    xydatah_co{1,k} = ref_co; 
    xydatah_co{2,k} = est_co; 

    [est,~,~] = Sinv(nvar, 5,i,j,2, 10);
    [~,col] = find(est==0); %finds zeros
    est1 = fliplr(est(:,1:col(1))); %set straight first part
    est2 = est(:,col(1):end);
    sz1 = size(est1)-1;sz2 = size(est2)-1;
    ref_est1 = (0:sz1(2))*delta_targetv_co;
    ref_est2 = (0:sz2(2))*delta_targetv_co;
    est_co = [est1(:)' est2(:)'];
    ref_est1 = repelem(ref_est1,sz1(1)+1,1);
    ref_est2 = repelem(ref_est2,sz1(1)+1,1);
    ref_co = [ref_est1(:)',ref_est2(:)'];
    xydatav_co{1,k} = ref_co; 
    xydatav_co{2,k} = est_co;  
end

%%
tpath = 37; %tamaño del torbulador 
zpath = @(z)[z;(tpath+z)]/L;
%% Produce a fit. Define parameters to optimize.
kp = optimvar('kp',1,'LowerBound',0,'UpperBound',100);
q = optimvar('q',1,'LowerBound',0,'UpperBound',100);
l0final = l0./kp; L0final = l0./q;
%% Data
sols = cell(1,M*M);
for k=1:100
    ydata = xydatah_co{2,k}';x = xydatah_co{1,k}';

    eval([ 'fun = @(kp,q) Stot(x,0,kp,q,1/3,1,' num2str(Fm_co) ',' num2str(deltamh_co) ...
    ',' num2str(Rmlab) ',.0001, [1 ; zpath(' num2str(47) ')]);']);
    % optimization problem
    response_co = fcn2optimexpr(fun,kp,q);
    obj = sum((response_co - ydata).^2);
    lsqproblem = optimproblem("Objective",obj);
    init0.kp = .8; init0.q = .1;
    %
    show(lsqproblem)
%
    tic,
    [sol,fval] = solve(lsqproblem,init0);
    toc
    sols{k} = sol; 
end
%%
kp = 8.4947;
q = 0.1866;

result = (20e-6/0.085)^2 * aoa / ...
    integralvar(2, Fm_co * kp^2, Rmlab * kp, 'divergent', 0.0001, 1/3, q, [1; zpath(47)])/ ...
    variance_constant(1/3, 0.84, Rmlab * kp, kmlab * kp, q);

save('C:\Yop\Seminario\Seminario_matlab_results\Optimizacion_mia\resultado_c^2n_V1725dT20_incoh.mat', 'result')