%             FEM CODE FOR LINEAR ELASTICITY 
% get wall time
time0 = datetime("now");
format long;

% ---open input & output file

global in;
in = fopen('mesh_1.inp', 'r');

global out;
out = fopen('result_1.out','w');

%--- Input FEM data:
[npoin,nelem,nvfix,ntype,nnode,ndofn,ndime,ngaus, ...
nstre,nmats,nprop,lnods,matno,coord,props,nofix, ...
iffix,fixed] = input_fem_elast( );
ntotv = npoin*ndofn;
[posgp, weigp]  = gauss(ngaus,nnode);

%--- Form global stiffnes matrix (lhs):
[gstif] = stiffness(npoin,nelem,nnode,nstre,ndime,ndofn, ...
    ngaus,ntype,lnods,matno,coord,props,posgp,weigp);

%--- Force vector & Boundary Conditions:
[gforce] = loads(npoin,nelem,ndofn,nnode,ngaus,ndime,posgp, ...
    weigp,lnods,coord);
[gstif,gforce]  = boundary_cond(npoin,nvfix,nofix,iffix, ...
    fixed,ndofn,gstif,gforce);

%--- Solve displacements :
asdis = gstif\gforce;

%--- calculate stresses:
[elem_stres] =  stress(asdis,nelem,npoin,nnode,ngaus,nstre,props, ...
    ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp);

%--- output results:
output(npoin,nelem,nnode,lnods,coord,ndofn,ngaus,nstre,...
    asdis,elem_stres);

%--- end of execution:
disp('done');
compute_time = datetime('now') - time0;

