D_2 = 0.03 
D_3 = 11 
nex_gef = 0.2 
nex_intr = 0.12 
hyd_42 = 1 
att_42bc = 0.266 
att_42 = 0.28 
ext_42 = 1 
att_b = 0.2667 
det_b = 0.35
att_24 = 0.00297 
det_24 = 0.35
R_ = 3.95 
N_42 = 3000
N_B = 6500
N_24 = 1000


function delta_mt = delta_mt(md, mt, mbg, thet, phi) 
delta_mt = (mbg*nex_gef + nex_intr)*md - hyd_42*mt + att_42bc*mbg*cd + D_2*laplacian(mt, [R_, thet, phi], %)

function delta_md = delta_md(mbg, md, mt, cd, thet, phi)
delta_md = -(mbg*nex_gef + nex_intr)*md + hyd_42*mt + att_42*cd - ext_42*md + D_2*laplacian(md, [R_, thet, phi], %)

function delta_mb = delta_mb(mt, mb, mb, cg, thet, phi)
delta_mb = att_b*mt*cb - det_b*mb*cg - att_24*mb*cb + det_24*mbg + D_2*laplacian(mb, [R_, thet, phi], %)

function delta_bg = delta_bg(mb, cg, mbg, thet, phi)
delta_bg = att_24*mb*cb - det_24 - det_24*mbg + D_2*laplacian(mbg, [R_, thet, phi], %)

function delta_cg = delta_cg(cd, r, thet, phi) 
delta_cg = D_3 * laplacian(cd, [r, thet, phi], %)

function delta_cb = delta_cb(cb, r, thet, phi) 
delta_cb = D_3 * laplacian(cb, [r, thet, phi], %)

function delta_cg = delta_cg(cg, r, thet, phi) 
delta_cg = D_3 * laplacian(cg, [r, thet, phi], %)

%boundary conditions 
function bound_cd = bound_cd(cd, mbg, md) 
bound_cd = (-(att_42bc * mbg + att_42) * cd + ext_42 * md) / D_3 

function bound_cb = bound_cb(mt, cb, mb) 
bound_cb = (-(att_b * mt* cb) + det_b * mb) / D_3 

function bound_cg = bound_cg(mb, cg, mbg) 
bound_cg = (-(att_24 * mb * cg) + det_24 * mbg) / D_3

