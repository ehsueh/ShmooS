% Create a 3D shperical space
% linalg::ogCoordTab[Spherical[RightHanded], Transformation](r, thet, phi);

% Define variables from Shooms Var File.
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
R_ = 3.95 ;
N_42 = 3000;
N_B = 6500;
N_24 = 1000;

% Define range
% Restrict the radius of shperical to be less than R_
r=0:0.1:R_;
thet=0:pi/4:2pi;
pi=-pi/2:pi/4:pi/2;

% concentration as a function of r, theta and phi

fucntion f_mt = f_mt(r, t, p)
f_mt = MT(r, t, p)

f_md = f_mt 
f_mb = f_mt
f_mbg = f_mt
f_cd = f_mt
f_cb = f_mt
f_cg = f_mt

% define matrices for mt, bd, mb, mbg, cd, cb, cg



% set timer



function delta_mt = delta_mt(md, mt, mbg, thet, phi) 
delta_mt = (mbg*nex_gef + nex_intr)*md - hyd_42*mt + att_42bc*mbg*cd + D_2*laplacian(f_mt, [R_, thet, phi], Spherical[RightHanded]);

function delta_md = delta_md(mbg, md, mt, cd, thet, phi)
delta_md = -(mbg*nex_gef + nex_intr)*md + hyd_42*mt + att_42*cd - ext_42*md + D_2*laplacian(f_md, [R_, thet, phi], Spherical[RightHanded]);

function delta_mb = delta_mb(mt, mb, mb, cg, thet, phi)
delta_mb = att_b*mt*cb - det_b*mb*cg - att_24*mb*cb + det_24*mbg + D_2*laplacian(f_mb, [R_, thet, phi], Spherical[RightHanded] );

function delta_bg = delta_bg(mb, cg, mbg, thet, phi)
delta_bg = att_24*mb*cb - det_24 - det_24*mbg + D_2*laplacian(f_mbg, [R_, thet, phi], Spherical[RightHanded]);

function delta_cg = delta_cg(cd, r, thet, phi) 
delta_cg = D_3 * laplacian(f_cd, [r, thet, phi], Spherical[RightHanded]);

function delta_cb = delta_cb(cb, r, thet, phi) 
delta_cb = D_3 * laplacian(f_cb, [r, thet, phi], Spherical[RightHanded]);

function delta_cg = delta_cg(cg, r, thet, phi) 
delta_cg = D_3 * laplacian(f_cg, [r, thet, phi], Spherical[RightHanded]);

%boundary conditions 
function bound_cd = bound_cd(cd, mbg, md) 
bound_cd = (-(att_42bc * mbg + att_42) * cd + ext_42 * md) / D_3 ;

function bound_cb = bound_cb(mt, cb, mb) 
bound_cb = (-(att_b * mt* cb) + det_b * mb) / D_3;

function bound_cg = bound_cg(mb, cg, mbg) 
bound_cg = (-(att_24 * mb * cg) + det_24 * mbg) / D_3;


% draw per time t

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

