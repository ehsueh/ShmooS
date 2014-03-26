% Create a 3D shperical space
% linalg::ogCoordTab[Spherical[RightHanded], Transformation](r, thet, phi);

% Define variables from Shooms Var File.

%% initial concentration
i_mt = 1;
i_md = 1;
i_mb = 1;
i_mbg = 1;
i_cd = 1;
i_cb = 1;
i_cg = 1;

%% time interval definition: 1 is the unit of each time step
% t - total time steps
time = 1;

%% display time interval
display_t = 1;


%% others
D_2 = 0.03; 
D_3 = 11 ;
nex_gef = 0.2 ;
nex_intr = 0.12; 
hyd_42 = 1 ;
att_42bc = 0.266; 
att_42 = 0.28; 
ext_42 = 1 ;
att_b = 0.2667 ;
det_b = 0.35;
att_24 = 0.00297; 
det_24 = 0.35;

N_42 = 3000;
N_B = 6500;
N_24 = 1000;

%% coordinates variables

% Radius variables: 
% i_r - inital value, R_ - max radius, r - interval, nr - number of r points
i_r = 0;
R_ = 3.95 ;
r_int = 0.01;
nr = (R_ - i_r) / r_int + 1;

% theta variables:
i_thet = 0;
THET_ = 2 * pi;
thet_int = pi / 4;
nthet = (THET_ - i_thet) / thet_int + 1;

% phi variables:
i_phi = - pi / 2;
PHI_ = pi / 2;
phi_int = pi / 4;
nphi = (PHI_ - i_phi) / phi_int + 1;

%% matrix definition and inquiry functions
% define matrices and functions for mt, bd, mb, mbg, cd, cb, cg
syms f r t p;

% initial materix to get the functions
cons = 0.01;
mat_mt = rand(nthet, nphi) * i_mt;
f_mt = t + p;

mat_md = rand(nthet, nphi) * i_md;
f_md = t + p;

mat_mb = ones(nthet, nphi) * i_mb;
f_mb = t + p;

mat_mbg = ones(nthet, nphi) * i_mbg;
f_mbg = t + p;
% f_mbg = @(t, p) mat_mbg((t - i_thet) / thet_int + 1, (p - i_phi) / phi_int + 1);

mat_cd = ones(nr, nthet, nphi) * i_cd;
f_cd = r + t + p;
% f_cd = @(r, t, p) mat_cd((r - i_r) / r_int + 1, (t - i_thet) / thet_int + 1, (p - i_phi) / phi_int + 1);

mat_cb = ones(nr, nthet, nphi) * i_cb;
f_cb = r + t + p;
% f_cb = @(r, t, p) mat_cb((r - i_r) / r_int + 1, (t - i_thet) / thet_int + 1, (p - i_phi) / phi_int + 1);

mat_cg = ones(nr, nthet, nphi) * i_cg;
f_cg = r + t + p;
% f_cg = @(r, t, p) mat_cg((r - i_r) / r_int + 1, (t - i_thet) / thet_int + 1, (p - i_phi) / phi_int + 1);

%% define the spherical laplacian

% spherical laplacian               

% spherical laplacian 

laplacian = 1/r^2 * diff((r^2 * diff(f, r)), r) + 1/r^2/sin(t) * diff(sin(t) * diff(f, t), t) + 1/r^2/sin(t)^2 * diff(f, p, 2);
laplacian2 = 1/r^2/sin(t) * diff(sin(t) * diff(f, t), t) + 1/r^2/sin(t)^2 * diff(f, p, 2);

%% delta functions: change in concentration per time unit

%%% those functions needs to be adjusted
% delta_mt = @(r, t, p)(f_mbg(r, t, p) * nex_gef + nex_intr) * f_md(r, t, p) - hyd_42 * f_mt(r, t, p) + att_42bc * f_mbg(r, t, p) * f_cd(r, t, p) + D_2 * laplacian2(f_mt, R_, t, p);
% 
% delta_md = @(r, t, p) -(f_mbg(r, t, p) * nex_gef + nex_intr) * f_md(r, t, p) + hyd_42 * f_mt(r, t, p) + att_42 * f_cd(r, t, p) - ext_42 * f_md(r, t, p) + D_2 * laplacian2(f_md, R_, t, p);
% 
% delta_mb = @(r, t, p) att_b * f_mt(r, t, p) * f_cb(r, t, p) - det_b * f_mb(r, t, p) * f_cg(r, t, p) - att_24 * f_mb(r, t, p) * f_cb(r, t, p) + det_24 * f_mbg(r, t, p) + D_2 * laplacian2(f_mb, R_, t, p);
% 
% delta_bg = @(r, t, p) att_24 * f_mb(r, t, p) * f_cb(r, t, p) - det_24 - det_24 * f_mbg(r, t, p) + D_2 * laplacian2(f_mbg, R_, t, p);

delta_mt = (f_mbg(t, p) * nex_gef + nex_intr) * f_md(t, p) - hyd_42 * f_mt(t, p) + att_42bc * f_mbg(t, p) * subs(f_cd, r, R_) + D_2 * subs(laplacian2, [f r], [f_mt R_]);

delta_md = -(f_mbg(t, p) * nex_gef + nex_intr) * f_md(t, p) + hyd_42 * f_mt(t, p) + att_42 * subs(f_cd, r, R_) - ext_42 * f_md(t, p) + D_2 * subs(laplacian2, [f r], [f_md R_]);

delta_mb =  att_b * f_mt(t, p) * subs(f_cb, r, R_) - det_b * f_mb(t, p) * subs(f_cg, r, R_) - att_24 * f_mb(t, p) * subs(f_cb, r, R_) + det_24 * f_mbg(t, p) + D_2 * subs(laplacian2, [f r], [f_mb R_]);

delta_mbg = att_24 * f_mb(t, p) * subs(f_cb, r, R_) - det_24 - det_24 * f_mbg(t, p) + D_2 * subs(laplacian2, [f r], [f_mbg R_]);

delta_cd = D_3 * laplacian(f_cd, r, t, p);

delta_cb = D_3 * laplacian(f_cb, r, t, p);

delta_cg = D_3 * laplacian(f_cg, r, t, p);

%% boundary conditions 
bound_cd = @(cd, mbg, md) (-(att_42bc * mbg + att_42) * cd + ext_42 * md) / D_3 ;

bound_cb = @(mt, cb, mb) (-(att_b * mt* cb) + det_b * mb) / D_3;

bound_cg = @(mb, cg, mbg) (-(att_24 * mb * cg) + det_24 * mbg) / D_3;

%% update functions


%% now start to run for time perid t
for tstep = 0 : time - 1
    
    % plot if it is at the right interval
%     if not mod(tstep, t_interval)
%     end 
    
    % creat an update materix and added to the existing one at the end of
    % interval
    % inital the delta matrices
%     mat_delta_mt = zeros(nr, nthet, nphi);
%     mat_delta_md = zeros(nr, nthet, nphi);
%     mat_delta_mb = zeros(nr, nthet, nphi);
%     mat_delta_mbg = zeros(nr, nthet, nphi);
%     mat_delta_cd = zeros(nr, nthet, nphi);
%     mat_delta_cb = zeros(nr, nthet, nphi);
%     mat_delta_cg = zeros(nr, nthet, nphi);
  
    % at 0 step, delta 1 and function are already there. 
    % calculate the function1
%     f_mt = @(r, t, p) delta_mt(R_, t, p) + f_mt(R_, t, p);
%     f_md = @(r, t, p) delta_md(R_, t, p) + f_md(R_, t, p);
%     f_mb = @(r, t, p) delta_mb(R_, t, p) + f_mb(R_, t, p);
%     f_mbg = @(r, t, p) delta_mbg(R_, t, p) + f_mbg(R_, t, p);
    f_mt = delta_mt(t, p) + f_mt(t, p);
    f_md = delta_md(t, p) + f_md(t, p);
    f_mb = delta_mb(t, p) + f_mb(t, p);
    f_mbg = delta_mbg(t, p) + f_mbg(t, p);
    f_cd = delta_cd(r, t, p) + f_cd(r, t, p);
    f_cb = delta_cb(r, t, p) + f_cb(r, t, p);
    f_cg = delta_cg(r, t, p) + f_cg(r, t, p);
    
    delta_mt = (f_mbg(t, p) * nex_gef + nex_intr) * f_md(t, p) - hyd_42 * f_mt(t, p) + att_42bc * f_mbg(t, p) * subs(f_cd, r, R_) + D_2 * subs(laplacian2, [f r], [f_mt R_]);

    delta_md = -(f_mbg(t, p) * nex_gef + nex_intr) * f_md(t, p) + hyd_42 * f_mt(t, p) + att_42 * subs(f_cd, r, R_) - ext_42 * f_md(t, p) + D_2 * subs(laplacian2, [f r], [f_md R_]);

    delta_mb =  att_b * f_mt(t, p) * subs(f_cb, r, R_) - det_b * f_mb(t, p) * subs(f_cg, r, R_) - att_24 * f_mb(t, p) * subs(f_cb, r, R_) + det_24 * f_mbg(t, p) + D_2 * subs(laplacian2, [f r], [f_mb R_]);

    delta_mbg = att_24 * f_mb(t, p) * subs(f_cb, r, R_) - det_24 - det_24 * f_mbg(t, p) + D_2 * subs(laplacian2, [f r], [f_mbg R_]);

    delta_cd = D_3 * laplacian(f_cd, r, t, p);

    delta_cb = D_3 * laplacian(f_cb, r, t, p);

    delta_cg = D_3 * laplacian(f_cg, r, t, p);

end

subs(f_mt, [t, p], [pi/4 pi/4])
 
% syms r t p;
% mt = r*2 + t^2 + p;
% 
% % create vectors
% thet = linspace(0, 2*pi);
% phi = linspace(-pi/2, pi/2);
% 
% % create meshgrid for inputs
% [thet, phi] = meshgrid(thet, phi);
% 
% % define radius
% r = R_;
% 
% % convert to cartesian coordinate
% [X, Y, Z] = sph2cart(thet, phi, r);
% for i=1:length(X)
% 	[THETA, PHI, R] = cart2sph(X(i), Y(i), Z(i)); 
% 	C(i) = subs(mt, [r t p], [R(1) THETA(1) PHI(1)]);
% end
% 
% % plot
% scatter3(X, Y, Z, 1, C);


