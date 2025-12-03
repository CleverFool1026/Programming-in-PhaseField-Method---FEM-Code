function[elem_stres] =  stress(asdis,nelem,npoin,nnode,ngaus,...
    nstre,props,ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp)
format long;

%--- Number of integration

ngaus2 = ngaus;
if(nnode == 3)
    ngaus2 = 1;
end
mgaus = ngaus*ngaus2;

for ielem = 1:nelem
%--- Material Parameters and Elasticity Matrix:

    mtype = matno(ielem);
    dmatx = modps(mtype,ntype,nstre,props);
    poiss = props(mtype,2);
    
    %--- Nodal Displacements:
    
    for inode = 1:nnode
        lnode = lnods(ielem,inode);
        for idofn = 1:ndofn
    
            nposn = (lnode-1)*ndofn+idofn;
            eldis(idofn,inode) = asdis(nposn);
            elcod(idofn,inode) = coord(lnode,idofn);
        end
    end
    
    %--- Integrate Stresses:
    
    kgasp = 0;
    for igaus = 1:ngaus
        for jgaus = 1:ngaus2
    
            kgasp = kgasp+1;
            exisp = posgp(igaus);
            etasp = posgp(jgaus);
    
            if(ngaus2 == 1)
                etasp = posgp(ngaus+igaus)
            end
    
            [shape,deriv] = sfr2(exisp,etasp,nnode);
            [cartd,djacb,gpcod] = jacob2(ielem,elcod,kgasp,shape,deriv,nnode,ndime);
            [bmatx] = bmats(cartd,shape,inode);
    
            %--- Calculate the strains:
    
            for istre = 1:nstre
                stran(istre) =0.0;
                for inode = 1:nnode
                    for idofn = 1:ndofn
                        ievab = (inode-1)*ndofn+idofn;
                        stran(istre) = stran(istre)+bmatx(istre,ievab)*eldis(idofn,inode);
                    end
                end
            end

%--- Calculate stresses:

            for istre = 1:nstre
                stres(istre) = 0.0;
                for jstre = 1:nstre
                    stres(istre) = stres(istre)+dmatx(istre,jstre)*stran(jstre);
                end
            end
            
            if(ntype == 1)
                stres(4) = 0.0;
            end
            if(ntype == 2)
                stres(4) = poiss*(stres(1)+stres(2));
            end
            
            for istre = 1:nstre+1
                elem_stres(ielem,kgasp,istre) = stres(istre);
            end


        end % igaus
    end % jgaus

end % ielem

end