function [gstif] = stiffness(npoin,nelem,nnode,nstre,...
    ndime,ndofn,ngaus,ntype,lnods,matno,coord,...
    props,posgp,weigp)
format long;

%--- Initialize the global stiffness:
nevab = nnode * ndofn;
ngaus2 = ngaus;
if(nnode == 3)
    ngaus2 = 1;
end

ntotv = npoin * ndofn;
gstif = zeros(ntotv, ntotv);

% ---Element stifness & loads:
for ielem = 1:nelem %大循环

%--- Initialize element stiffness:
for ievab = 1:nevab
    for jevab = 1:nevab
        estif(ievab, jevab) = 0.0;
    end
end

%--- Form elasticity matrix:
mtype = matno(ielem);
dmatx =modps(mtype,ntype,nstre,props);

%--- Coordinates of element nodes:
for inode = 1:nnode
    lnode = lnods(ielem,inode);
    for idime = 1:ndime
        elcod(idime,inode) = coord(lnode,idime);
    end
end

%--- Integrate the element stiffness :
kgasp = 0;
for igaus = 1:ngaus
    exisp = posgp(igaus);
    for jgaus = 1:ngaus2
        etasp  = posgp(jgaus);
        if(nnode == 3)
            etasp = posgp(ngaus + igaus);
        end

        kgasp = kgasp + 1;

        [shape, deriv] = sfr2(exisp, etasp, nnode);
        [cartd,djacb,gpcod] = jacob2(ielem,elcod,kgasp, ...
            shape,deriv,nnode,ndime);
        [bmatx] = bmats(cartd, shape, inode);
        [dbmat] = dbe(nevab, nstre, bmatx, dmatx);

        dvolu = djacb*weigp(igaus)*weigp(jgaus);

        if(node == 3)
            dvolu = djacb*weigp(igaus);
        end

%---Form element stiffness:
        for ievab = 1:nevab
            for jevab = 1:nevab
                for istre = 1:nstre
                    estif(ievab,jevab) = estif(ievab,jevab)+bmatx(istre,ievab)*...
                    dbmat(istre,jevab)*dvolu;
                end
            end
        end

end%jgaus
end%igaus

%--- Form Global stiffness matrix:
for inode = 1:nnode
    lnode = lnods(ielem,inode);
    for idofn = 1:ndofn
        itotv = (lnode-1)*ndofn+idofn;
        ievab = (inode-1)*ndofn+idofn;
        for jnode = 1:nnode
            knode = lnods(ielem,jnode);
            for jdofn = 1:ndofn
                jtotv = (knode-1)*ndofn+jdofn;
                jevab = (jnode-1)*ndofn+jdofn;
                
                gstif(itotv,jtotv) = gstif(itotv,jtotv)+...
                estif(ievab,jevab);
            end
        end
    end
end

end%ielem
end