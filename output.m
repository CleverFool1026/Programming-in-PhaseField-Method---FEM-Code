function[ ] = output(npoin,nelem,nnode,lnods,coord,ndofn,...
    ngaus,nstre,asdis,elem_stres)
format long;

global out;

fprintf(out,'\n');
fprintf(out,'********************\n');
fprintf(out,'* OUTPUTS *\n');
fprintf(out,'********************\n');

%--- Print Nodal Displacements:
fprintf(out,'\n');
fprintf(out,'Displacements:\n');
fprintf(out,'Node Number X-Disp Y-Disp\n');

for ipoin = 1:npoin
    fprintf(out,'%5d',ipoin);
    for idofn = 1:ndofn
        itotv = (ipoin-1)*ndofn+idofn;
        fprintf(out,'%14.6e',asdis(itotv));
    end
    fprintf(out,'\n');
end

%--- Print Stress values:

fprintf(out,'\n');
fprintf(out,'Stress at elements\n');
fprintf(out,'gauss-point sigma-11sigma-sigma-12 sigma-33\n');

%--- Number of integration

ngaus2 = ngaus;
if(nnode == 3)
    ngaus2 = 1;
end


for ielem = 1:nelem
    fprintf(out,'\n');
    fprintf(out,'Element No: %5d\n',ielem);

    kgasp = 0;
    for igaus = 1:ngaus
        for jgaus = 1:ngaus2
            kgasp = kgasp +1;

            fprintf(out,'%d %14.6e %14.6e %14.6e%14.6e\n',kgasp,...
                elem_stres(ielem,kgasp,1),elem_stres(ielem,kgasp,2),...
                elem_stres(ielem,kgasp,3),elem_stres(ielem,kgasp,4));

        end
    end
end % ielem

%--- prepare results for graphical output:
%-- deformed mesh:
facto = 3.0;

for ipoin = 1:npoin
    for idofn  = 1:ndofn
        itotv = (ipoin-1)*ndofn+idofn;
        disp_cord(ipoin,idofn) = coord(ipoin,idofn) +...
            facto * asdis(itotv);
    end
end

%--- Extraplot stresses from integrationpoints
%--- and average over entire mesh
%--number of connection of nodes:
for ipoin = 1:npoin
    node_con(ipoin) = 0;
    for ielem = 1:nelem
        for inode = 1:nnode
            lnode = lnods(ielem,inode);
            if(lnode == ipoin)
                node_con(ipoin) = node_con(ipoin) +1;
            end
        end
    end
end %ipoin

%-- initialize nodal_stress

for ipoin = 1:npoin
    for istre = 1:nstre
        node_stres(ipoin,istre) = 0.0;
    end
end

for ielem = 1:nelem

    for istre = 1:nstre
        ave_stres(istre) = 0.0;

    end

    kgasp = 0.0;
    for igaus = 1:ngaus
        for jgaus = 1:ngaus2
            kgasp = kgasp +1;

            for istre = 1:nstre
                ave_stres(istre) = ave_stres(istre)+elem_stres(ielem,kgasp,istre);
            end
        end
    end

    for inode = 1:nnode
        lnode = lnods(ielem,inode);
        for istre = 1:nstre
            node_stres(lnode,istre) = node_stres(lnode,istre) + ave_stres(istre)/kgasp;
        end
    end
end %ielem

for ipoin = 1:npoin
    for istre = 1:nstre
        node_stres(ipoin,istre) = node_stres(ipoin,istre)/node_con(ipoin);
    end
end

%-- switch order of element connectivity if nnode = 8
if(nnode == 8 )
    for ielem = 1:nelem
        for inode = 1:nnode
            dummy(inode) = lnods(ielem,inode);
        end
        lnods(ielem,1) = dummy(1);
        lnods(ielem,2) = dummy(3);
        lnods(ielem,3) = dummy(5);
        lnods(ielem,4) = dummy(7);
        lnods(ielem,5) = dummy(2);
        lnods(ielem,6) = dummy(4);
        lnods(ielem,7) = dummy(6);
        lnods(ielem,8) = dummy(8);
    end
end %if

%--- output to vtk file.
istep = 1;

write_vtk_fem(npoin,nelem,nnode,lnods,...
    coord,istep,node_stres)
end 