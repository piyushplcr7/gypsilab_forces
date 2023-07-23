function found_indices = findIndicesInBigVec(i,j,I0,I1,I2,I3,J0,J1,J2,J3,Nqud)
    % inflated indices vectors
    inflatedI = [repelem(I0,Nqud(1),1);...
                 repelem(I1,Nqud(2),1);...
                 repelem(I2,Nqud(3),1);...
                 repelem(I3,Nqud(4),1)];

    inflatedJ = [repelem(J0,Nqud(1),1);...
                 repelem(J1,Nqud(2),1);...
                 repelem(J2,Nqud(3),1);...
                 repelem(J3,Nqud(4),1)];

    fi = find(inflatedI == i);
    fj = find(inflatedJ(fi) == j);

    found_indices = fi(fj);

    % Checking if the found indices are right
    checkI = inflatedI(found_indices);
    checkJ = inflatedJ(found_indices);
    assert(size(unique(checkI),1) == 1);
    assert(size(unique(checkJ),1) == 1);

    assert(checkI(1) == i);
    assert(checkJ(1) == j);

end