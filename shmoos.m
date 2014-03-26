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
t = 100;

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
r_int = 0.1;
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
mat_mt = ones(nr, nthet, nphi) * i_mt;
f_mt = @(r, t, p) mat_mt(r, t, p);

mat_md = ones(nr, nthet, nphi) * i_md;
f_md = @(r, t, p) mat_md(r, t, p);

mat_mb = ones(nr, nthet, nphi) * i_mb;
f_mb = @(r, t, p) mat_mb(r, t, p);

mat_mbg = ones(nr, nthet, nphi) * i_mbg;
f_mbg = @(r, t, p) mat_mbg(r, t, p);

mat_cd = ones(nr, nthet, nphi) * i_cd;
f_cd = @(r, t, p) mat_cd(r, t, p);

mat_cb = ones(nr, nthet, nphi) * i_cb;
f_cb = @(r, t, p) mat_cb(r, t, p);

mat_cg = ones(nr, nthet, nphi) * i_cg;
f_cg = @(r, t, p) mat_cg(r, t, p);

%% delta functions: change in concentration per time unit

%%% those functions needs to be adjusted
delta_mt = @(md, mt, mbg, thet, phi)(mbg*nex_gef + nex_intr)*md - hyd_42*mt + att_42bc*mbg*cd + D_2*laplacian(f_mt, [R_, thet, phi], Spherical[RightHanded]);

delta_md = @(mbg, md, mt, cd, thet, phi) -(mbg*nex_gef + nex_intr)*md + hyd_42*mt + att_42*cd - ext_42*md + D_2*laplacian(f_md, [R_, thet, phi], Spherical[RightHanded]);

delta_mb = @(mt, mb, mb, cg, thet, phi) att_b*mt*cb - det_b*mb*cg - att_24*mb*cb + det_24*mbg + D_2*laplacian(f_mb, [R_, thet, phi], Spherical[RightHanded] );

delta_bg = @(mb, cg, mbg, thet, phi) att_24*mb*cb - det_24 - det_24*mbg + D_2*laplacian(f_mbg, [R_, thet, phi], Spherical[RightHanded]);

delta_cg = @(cd, r, thet, phi) D_3 * laplacian(f_cd, [r, thet, phi], Spherical[RightHanded]);

delta_cb = @(cb, r, thet, phi) D_3 * laplacian(f_cb, [r, thet, phi], Spherical[RightHanded]);

delta_cg = @(cg, r, thet, phi) D_3 * laplacian(f_cg, [r, thet, phi], Spherical[RightHanded]);

%% boundary conditions 
bound_cd = @(cd, mbg, md) (-(att_42bc * mbg + att_42) * cd + ext_42 * md) / D_3 ;

bound_cb = @(mt, cb, mb) (-(att_b * mt* cb) + det_b * mb) / D_3;

bound_cg = @(mb, cg, mbg) (-(att_24 * mb * cg) + det_24 * mbg) / D_3;

%% update functions
   

% create vectors 
thet = linspace(0, 2*pi);
phi = linspace(-pi/2, pi/2);

% create meshgrid for inputs
[thet, phi] = meshgrid(thet, phi);

% define radius
r = 1

% convert to cartesian coordinate
[x, y, z] = sph2cart(thet, phi, r);

% plot
plot(x, y, z);

