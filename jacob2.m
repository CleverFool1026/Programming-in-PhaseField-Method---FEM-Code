function[cartd, djacb, gpcod] = jacob2(ielem,elcod,kgasp,shape,...
deriv,nnode,ndime)
format long;

%--- Guass point coordinates:
for idime = 1:ndime
    gpcod(idime, kgasp) = 0.0;
    for inode = 1:nnode
        gpcod(idime,kgasp) = gpcod(idime, kgasp)+...
         elcod(idime,inode)*shape(inode);

    end
end

% --- Jacobian
for idime = 1:ndime
    for jdime = 1:ndime
        xjacm(idime,jdime) = 0.0;
        for inode  =1:nnode
            xjacm(idime,jdime) = xjacm(idime,jdime)+...
             deriv(idime,inode)*elcod(jdime,inode);
        end
    end
end

djacb = xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1);
if(djacb <= 0.0)
    disp('Program Terminated')
    disp('Zero or negative area for element:')
    disp(ielem);
    error('Program terminated zero or negative area')
end

%--- Cartesian derivatives:
xjaci(1,1) = xjacm(2,2)/djacb;
xjaci(2,2) = xjacm(1,1)/djacb;
xjaci(1,2) = -xjacm(1,2)/djacb;
xjaci(2,1) = -xjacm(2,1)/djacb;

for idime = 1:ndime
    for inode = 1:nnode
        cartd(idime,inode) = 0.0;
        for jdime = 1:ndime
            cartd(idime,inode) = cartd(idime,inode)+ ... 
                xjaci(idime,jdime)*deriv(jdime,inode);
        end
    end
end
end