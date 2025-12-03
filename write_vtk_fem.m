function [] = write_vtk_fem(npoin, nelem, nnode, lnods, coord, istep, cont1)

format long;

%-- open file
fname=sprintf('time_%d.vtk', istep);
out=fopen(fname,'w');

%-- start writing:

% header
fprintf(out,'# vtk DataFile Version 2.0\n');
fprintf(out,'time_10.vtk\n');
fprintf(out,'ASCII\n');
fprintf(out,'DATASET UNSTRUCTURED_GRID\n');

%write nodal coordinates:
fprintf(out,'POINTS %d float\n',npoin);
dummy=0.0;
for ipoin=1:npoin
    fprintf(out,'%14.6f %14.6f %14.6f\n', coord(ipoin,1),coord(ipoin,2), dummy);
end

%-- write element connectivity:
iconst1=nelem*(nnode+1);
fprintf(out,'CELLS %d %d\n',nelem, iconst1);
for ielem=1:nelem
    fprintf(out,'%d ',nnode);
    for inode=1:nnode
        fprintf(out,'%5d',(lnods(ielem, inode)-1));
    end
    fprintf(out,'\n');
end

%-- write cell types:
if(nnode == 8)
    ntype=23;
end
if(nnode == 4)
    ntype=9;
end
if(nnode == 3)
    ntype=5;
end

fprintf(out,'CELL_TYPES %d\n',nelem);

for i=1:nelem
    fprintf(out,'%2d\n', ntype);
end

%-- write Nodal scalar & vector values:
fprintf(out,'POINT_DATA %d\n',npoin);

%-- write stress values as scalar:
fprintf(out,'SCALARS sigma_xx float 1\n');
fprintf(out,'LOOKUP_TABLE default\n');

for ipoin=1:npoin
    fprintf(out,'%14.6e\n',cont1(ipoin,1));
end

fprintf(out,'SCALARS sigma_yy float 1\n');
fprintf(out,'LOOKUP_TABLE default\n');

for ipoin=1:npoin
    fprintf(out,'%14.6e\n',cont1(ipoin,2));
end

fprintf(out,'SCALARS sigma_xy float 1\n');
fprintf(out, 'LOOKUP_TABLE default\n');

for ipoin = 1:npoin
    fprintf(out,'%14.6e\n',cont1(ipoin,3));
end

fclose(out);

end %endfunction