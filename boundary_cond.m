function[gstif,gforce] = boundary_cond(npoin,nvfix,nofix,iffix,fixed,...
    ndofn,gstif,gforce)
format long;

ntotv = npoin*ndofn;

for ivfix = 1:nvfix
    lnode = nofix(ivfix);

    for idofn = 1:ndofn
        if(iffix(ivfix,idofn) == 1)
            itotv = (lnode - 1)*ndofn +idofn;

            for jtotv = 1:ntotv
                gstif(itotv,jtotv) = 0.0;
            end

            gstif(itotv,itotv) = 1.0;
            gforce(itotv) = fixed(ivfix,idofn);

        end
    end
end

end