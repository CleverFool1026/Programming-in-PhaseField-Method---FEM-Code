function[gforce] = loads(npoin,nelem,ndofn,nnode,ngaus,ndime, ...
    posgp,weigp,lnods,coord)
format long;

global in;
global out;

%--- Initilize global force vector (rhs) & element loads:
nevab = nnode * ndofn;
ntotv = npoin * ndofn;

for itotv = 1:ntotv
    gforce = zeros(itotv, 1);
end

for ielem = 1:nelem
    for ievab = 1:nevab
        eload(ielem,ievab) = 0.0;
    end
end

%--- Loading types:
dummy = fscanf(in, '%d %d', [2, 1]);
iplod = dummy(1);
nedge = dummy(2);

%--- Point forces:
if(iplod ~= 0)
    nplod = fscanf(in,'%d',[1,1]);
    for jplod = 1:nplod
        lodpt = fscanf(in,'%d',[1,1]);
        dummy = fscanf(in,'%lf %lf',[2,1]);
        
        for idofn = 1:ndofn
            point(idofn) = dummy(idofn);
        end
        
        for ielem = 1:nelem
            for inode = 1:nnode
                nloca = lnods(ielem,inode);
                if(lodpt == nloca)
                    for idofn = 1:ndofn
                        nposi = (inode-1)*ndofn+idofn;
                        eload(ielem,nposi) = point(idofn);
                    end
                end
            end
        end
    end    
end % if_iplod

%--- Distributed Forces:
if(nedge ~= 0)
    fprintf(out, 'Number of loaded edges:%d\n', nedge);
    fprintf(out, 'List of loaded edges and applied loads:\n');

    for iedge = 1:nedge
        nodeg = 3;
        if(nnode ~= 8)
            nodeg = 2;
        end

%--- Read loads:
        neass = fscanf(in,'%d',[1,1]);
        dummy = fscanf(in,'%d %d %d',[nodeg,1]);
        for iodeg = 1:nodeg
            noprs(iodeg) = dummy(iodeg);
        end
        
        for iodeg = 1:nodeg
            dummy = fscanf(in,'%lf %lf',[2,1]);
            for idofn = 1:ndofn
                press(iodeg,idofn) = dummy(idofn);
            end
        end

%--- Print:
        fprintf(out,'\n');
        fprintf(out,'%5d',neass);
        for iodeg = 1:nodeg
            fprintf(out,'%5d',noprs(iodeg));
        end
        fprintf(out,'\n');
        for iodeg = 1:nodeg
            for idofn = 1:ndofn
                fprintf(out,'%14.6e',press(iodeg,idofn));
            end
        end
        fprintf(out,'\n');

%-- End of reading
%--- Integrate along the edges:
        etasp = -1.0;

        for iodeg = 1:nodeg
            lnode = noprs(iodeg);
            for idime = 1:ndime
                elcod(idime,iodeg) = coord(lnode,idime);
            end
        end

        for igaus = 1:ngaus
            exisp = posgp(igaus);
            [shape, deriv] = sfr2(exisp,etasp,nnode);

            for idofn = 1:ndofn
                pgash(idofn) = 0.0;
                dgash(idofn) = 0.0;
                for iodeg = 1:nodeg
                    pgash(idofn) = pgash(idofn)+press(iodeg,idofn)*shape(iodeg);
                    dgash(idofn) = dgash(idofn)+elcod(idofn,iodeg)*deriv(1,iodeg);
                end
            end

            dvolu = weigp(igaus);
            pxcom = dgash(1)*pgash(2)-dgash(2)*pgash(1);
            pycom = dgash(1)*pgash(1)+dgash(2)*pgash(2);

            for inode = 1:nnode
                nloca = lnods(neass,inode);
                if(nloca == noprs(1))
                    jnode = inode+nodeg-1;
                    kount = 0;
                    for knode = inode:jnode
                        kount = kount+1;
                        ngash = (knode-1)*ndofn+1;
                        mgash = (knode-1)*ndofn+2;
                        if(knode> nnode)
                            ngash = 1;
                            mgash = 2;
                        end
                        eload(neass,ngash) = eload(neass,ngash)+pxcom*dvolu*shape(kount);
                        eload(neass,mgash) = eload(neass,mgash)+pycom*dvolu*shape(kount);
                    end
                end
            end

end %igauss
end % iedge

end%nedge

%--- Print Nodal Forces:
fprintf(out,'\n');
fprintf(out,'Nodal forces for elements:\n');
for ielem = 1:nelem

    fprintf(out, '\n');
    fprintf(out,'Element No: %d\n',ielem);
    for ievab = 1:nevab
        fprintf(out,'%14.6e',eload(ielem,ievab));
        if((nnode == 8) && (ievab == nevab /2))
            fprintf(out, '\n');
        end
    end
    fprintf(out,'\n');
end

%--- Generate global force vector:

for ielem = 1:nelem
    for inode = 1:nnode
        lnode = lnods(ielem,inode);
        for idofn = 1:ndofn
            itotv = (lnode-1)*ndofn+idofn;
            ievab = (inode-1)*ndofn+idofn;
            gforce(itotv) = gforce(itotv)+eload(ielem,ievab);
        end
    end
end

end