function [stress, nurbsStr, mises] = getStressBoundaryVisco(nurbsStr, ctrlDisp, ctrlTrac, deltaT, eTau, tau, visG, visK, collocNew)
% Calculates the surface stress and strain using the boundary element method.
%
% Inputs:
%   Str: Structure array containing boundary element data.
%   ctrlDisp: Displacement values at control points (n*3 matrix).
%   ctrlTrac: Traction values at control points (n*3 matrix).
%
% Outputs:
%   stress: Calculated surface stress.

numCollocT = size(ctrlTrac, 1);
stress = zeros(numCollocT, 6);
stressColl = zeros(size(collocNew.collocPts, 1), 6);
mises = zeros(size(collocNew.collocPts, 1), 1);
temp = 0;

E = nurbsStr(1).E;
nu = nurbsStr(1).nu;
G = E / 2 / (1 + nu);
K = E / (3 * (1 - 2 * nu));

for i = 1:length(nurbsStr)
    numColloc = nurbsStr(i).numCollocT;
    ind = temp + 1 : temp + numColloc;
    xi = nurbsStr(i).collocXi(:, 1);
    eta = nurbsStr(i).collocXi(:, 2);
    eBar = zeros(numColloc, 4);
    stressSub = zeros(numColloc, 6);
    misesSub = zeros(numColloc, 1);

    e0_bar1 = sum(nurbsStr(i).eBar1(:, 1 : 3), 2);
    e0_bar2 = sum(nurbsStr(i).eBar2(:, 1 : 3), 2);
    e0_bar1 = repmat(e0_bar1, [1, length(visG)]);
    e0_bar2 = repmat(e0_bar2, [1, length(visG)]);
    qK_i = eTau' .* nurbsStr(i).qK + K * visK' .* ((1 - eTau)' .* e0_bar1 ...
        + (1 - tau / deltaT .* (1 - eTau))' .* (e0_bar1 - e0_bar2));

    % attention: Str(i).eBar1(:, 1) - 1 / 3 * e0_bar1
    qG11_i = eTau' .* nurbsStr(i).qG11 + G * visG' .* ((1 - eTau)' .* ( nurbsStr(i).eBar1(:, 1) - 1 / 3 * e0_bar1) ...
        + (1 - tau / deltaT .* (1 - eTau))' .* (nurbsStr(i).eBar1(:, 1) - 1 / 3 * e0_bar1 - nurbsStr(i).eBar2(:, 1) + 1 / 3 * e0_bar2));

    qG22_i = eTau' .* nurbsStr(i).qG22 + G * visG' .* ((1 - eTau)' .* ( nurbsStr(i).eBar1(:, 2) - 1 / 3 * e0_bar1) ...
        + (1 - tau / deltaT .* (1 - eTau))' .* (nurbsStr(i).eBar1(:, 2) - 1 / 3 * e0_bar1 - nurbsStr(i).eBar2(:, 2) + 1 / 3 * e0_bar2));

    qG33_i = eTau' .* nurbsStr(i).qG33 + G * visG' .* ((1 - eTau)' .* ( nurbsStr(i).eBar1(:, 3) - 1 / 3 * e0_bar1) ...
        + (1 - tau / deltaT .* (1 - eTau))' .* (nurbsStr(i).eBar1(:, 3) - 1 / 3 * e0_bar1 - nurbsStr(i).eBar2(:, 3) + 1 / 3 * e0_bar2));

    % attention: / 2
    qG12_i = eTau' .* nurbsStr(i).qG12 + G * visG' .* ((1-eTau)' .* ( nurbsStr(i).eBar1(:, 4) / 2) ...
        + (1 - tau / deltaT .* (1 - eTau))' .* (nurbsStr(i).eBar1(:,4) - nurbsStr(i).eBar2(:,4)) / 2 );

    if nurbsStr(i).bouInd == 1
        for j = 1:numColloc
            qK = sum(qK_i(j, :), 2);
            qG11 = sum(qG11_i(j, :), 2);
            qG22 = sum(qG22_i(j, :), 2);
            qG33 = sum(qG33_i(j, :), 2);
            qG12 = sum(qG12_i(j, :), 2);
            [stressSub(j, :), eBar(j, :), misesSub(j)] = findBoundaryValue(ctrlDisp(ind,:), ctrlTrac(ind,:), ...
                [xi(j), eta(j)], nurbsStr(i), qK, qG11, qG22, qG33, qG12);
        end
    else
        for j = 1:numColloc
            qK = sum(qK_i(j, :), 2);
            qG11 = sum(qG11_i(j, :), 2);
            qG22 = sum(qG22_i(j, :), 2);
            qG33 = sum(qG33_i(j, :), 2);
            qG12 = sum(qG12_i(j, :), 2);
            [stressSub(j, :), eBar(j, :), misesSub(j)] = findBoundaryValue(ctrlDisp(ind,:), -ctrlTrac(ind,:), ...
                [xi(j), eta(j)], nurbsStr(i), qK, qG11, qG22, qG33, qG12);
        end
    end

    nurbsStr(i).eBar2 = nurbsStr(i).eBar1;
    nurbsStr(i).eBar1 = eBar;
    nurbsStr(i).qK = qK_i;
    nurbsStr(i).qG11 = qG11_i;
    nurbsStr(i).qG22 = qG22_i;
    nurbsStr(i).qG33 = qG33_i;
    nurbsStr(i).qG12 = qG12_i;

    temp = temp + numColloc;
    stress(ind, :) = stressSub;
    % mises(nurbsStr(i).globU, :) = misesSub + mises(nurbsStr(i).globU, :);
    stressColl(nurbsStr(i).globU, :) = stressSub + stressColl(nurbsStr(i).globU, :);
end

xx = [];
for i = 1 : length(nurbsStr)
    [xx] = [xx; nurbsStr(i).globU'];
end
yy = tabulate(xx(:));
zz = yy(:, 2);

for i = 1 : size(collocNew.collocPts, 1)
    stressColl(i,:) = stressColl(i,:) / zz(i);

    % mises(i,:)=mises(i,:)/zz(i);

    tempStress = [stressColl(i, 1),stressColl(i, 4),stressColl(i, 6);
        stressColl(i, 4),stressColl(i, 2),stressColl(i, 5);
        stressColl(i, 6),stressColl(i, 5),stressColl(i, 3)];
    sP = eig(tempStress);
    mises(i) = 1/sqrt(2)*sqrt((sP(1)-sP(2))^2 + (sP(2)-sP(3))^2 + (sP(3)-sP(1))^2);
end

end

function [stressOut, strain_bar, mises] = findBoundaryValue(ctrlDisp, ctrlTrac, paramCoords, nurbsStr, qK, qG11, qG22, qG33, qG12)
% Calculate the strain and stress at a given collocation point using the boundary values.
%
% Inputs:
%   ctrlDisp: Displacement values at control points (n*3 matrix).
%   ctrlTrac: Traction values at control points (n*3 matrix).
%   paramCoords: Parametric coordinates [xi, eta] of the collocation point.
%   nurbsStr: Structure containing boundary element data for a single patch.
%
% Outputs:
%   stressOut: Calculated stress at the collocation point.

p = nurbsStr.p;
q = nurbsStr.q;
knotVecU = nurbsStr.knotU;
knotVecV = nurbsStr.knotV;
weights = nurbsStr.weights;
numElem = nurbsStr.numElem;
xiU = paramCoords(1);
xiV = paramCoords(2);
E = nurbsStr.E;
nu = nurbsStr.nu;
G = E / 2 / (1 + nu);
K = E / (3 * (1 -2 * nu));

for elem = 1:numElem
    range = nurbsStr.elemRange(elem, :);
    globBasisIdx = nurbsStr.globIdx(elem, :);
    uConn = nurbsStr.uConn(elem, :);
    tConn = nurbsStr.tConn(elem, :);
    elCoords = nurbsStr.ctrlPts(globBasisIdx, 1:3);

    if range(1) <= xiU && xiU <= range(2) && range(3) <= xiV && xiV <= range(4)
        [XiU, XiV] = getXiParamter(xiU, xiV, range);
        [N, dNu, dNv] = NURBS2DBasisDers([XiU, XiV], p, q, knotVecU, knotVecV, weights');

        dxdxi = [dNu; dNv] * elCoords;
        Jxi = norm(dxdxi(1, :));
        Jeta = norm(dxdxi(2, :));
        normalQ = cross(dxdxi(1, :), dxdxi(2, :)) / norm(cross(dxdxi(1, :), dxdxi(2, :)));
        V1 = dxdxi(1, :) / Jxi;
        V2 = cross(normalQ, V1) / norm(cross(normalQ, V1));
        cosTheta = dot(V1, dxdxi(2, :) / Jeta);
        sinTheta = sqrt(1 - cosTheta^2);

        d_xi_d_x_bar = 1 / Jxi;
        d_xi_d_y_bar = -cosTheta / (Jxi * sinTheta);
        d_eta_d_x_bar = 0;
        d_eta_d_y_bar = 1 / Jeta / sinTheta;

        Ux = ctrlDisp(uConn, 1);
        Uy = ctrlDisp(uConn, 2);
        Uz = ctrlDisp(uConn, 3);
        dudPara = [dNu; dNv] * [Ux, Uy, Uz];
        du_dxi = dudPara(1, :);
        du_deta = dudPara(2, :);

        Tx = ctrlTrac(tConn, 1);
        Ty = ctrlTrac(tConn, 2);
        Tz = ctrlTrac(tConn, 3);
        Tpoint = N * [Tx, Ty, Tz];

        LTmatrix = [V1; V2; normalQ];
        Txyz_bar = LTmatrix * Tpoint';
        sigma_xz_bar = Txyz_bar(1);
        sigma_yz_bar = Txyz_bar(2);
        sigma_z_bar = Txyz_bar(3);

        epsilon_x_bar = dot(du_dxi, V1) * d_xi_d_x_bar + dot(du_deta, V1) * d_eta_d_x_bar;
        epsilon_y_bar = dot(du_dxi, V2) * d_xi_d_y_bar + dot(du_deta, V2) * d_eta_d_y_bar;
        epsilon_xy_bar = dot(du_dxi, V1) * d_xi_d_y_bar + dot(du_deta, V1) * d_eta_d_y_bar + ...
            dot(du_dxi, V2) * d_xi_d_x_bar + dot(du_deta, V2) * d_eta_d_x_bar;


        % sigma_x_bar = E / (1 - nu^2) * (epsilon_x_bar + nu * epsilon_y_bar) + nu / (1 - nu) * sigma_z_bar;
        % sigma_y_bar = E / (1 - nu^2) * (epsilon_y_bar + nu * epsilon_x_bar) + nu / (1 - nu) * sigma_z_bar;
        % sigma_xy_bar = E / (2 * (1 + nu)) * epsilon_xy_bar;
        % Viscoelastic
        % ε33
        epsilon_z_bar = (sigma_z_bar + (2 * G / 3 - K) * (epsilon_x_bar + epsilon_y_bar) ...
            + 2 * qG33 + qK ) / (K + 4 / 3 * G);

        strain_bar = [epsilon_x_bar; epsilon_y_bar; epsilon_z_bar; epsilon_xy_bar];

        epsilon_0_bar = epsilon_x_bar + epsilon_y_bar + epsilon_z_bar;

        % σ11, σ22, σ12
        sigma_x_bar = 2 * G * epsilon_x_bar + (K - 2 * G / 3) * epsilon_0_bar - 2 * qG11 - qK;
        sigma_y_bar = 2 * G * epsilon_y_bar + (K - 2 * G / 3) * epsilon_0_bar - 2 * qG22 - qK;
        sigma_xy_bar = G * epsilon_xy_bar - 2 * qG12;

        sigma_bar = [sigma_x_bar; sigma_y_bar; sigma_z_bar; sigma_xy_bar; sigma_yz_bar; sigma_xz_bar];

        L11 = LTmatrix(1, 1); L12 = LTmatrix(1, 2); L13 = LTmatrix(1, 3);
        L21 = LTmatrix(2, 1); L22 = LTmatrix(2, 2); L23 = LTmatrix(2, 3);
        L31 = LTmatrix(3, 1); L32 = LTmatrix(3, 2); L33 = LTmatrix(3, 3);

        QT_sigma_mat2 = ...
            [ L11^2, L21^2, L31^2, 2*L11*L21, 2*L21*L31, 2*L31*L11;
            L12^2, L22^2, L32^2, 2*L12*L22, 2*L22*L32, 2*L32*L12;
            L13^2, L23^2, L33^2, 2*L13*L23, 2*L23*L33, 2*L33*L13;
            L11*L12, L21*L22, L31*L32, L11*L22+L21*L12, L21*L32+L31*L22, L31*L12+L11*L32;
            L12*L13, L22*L23, L32*L33, L12*L23+L22*L13, L22*L33+L32*L23, L32*L13+L12*L33;
            L13*L11, L23*L21, L33*L31, L13*L21+L23*L11, L23*L31+L33*L21, L33*L11+L13*L31; ];

        sigma = QT_sigma_mat2 * sigma_bar;

        %         QT_sigma_mat = [
        %             L11^2, L12^2, L13^2, 2 * L11 * L12, 2 * L12 * L13, 2 * L13 * L11;
        %             L21^2, L22^2, L23^2, 2 * L21 * L22, 2 * L22 * L23, 2 * L23 * L21;
        %             L31^2, L32^2, L33^2, 2 * L31 * L32, 2 * L32 * L33, 2 * L33 * L31;
        %             L11 * L21, L12 * L22, L13 * L23, L11 * L22 + L12 * L21, L12 * L23 + L13 * L22, L13 * L21 + L11 * L23;
        %             L21 * L31, L22 * L32, L23 * L33, L21 * L32 + L22 * L31, L22 * L33 + L23 * L32, L23 * L31 + L21 * L33;
        %             L31 * L11, L32 * L12, L33 * L13, L31 * L12 + L32 * L11, L32 * L13 + L33 * L12, L33 * L11 + L31 * L13;
        %             ];
        %
        %         sigma = QT_sigma_mat \ sigma_bar;

        stress = [
            sigma(1), sigma(4), sigma(6);
            sigma(4), sigma(2), sigma(5);
            sigma(6), sigma(5), sigma(3)
            ];
        sPri = eig(stress);
        mises = 1/sqrt(2)*sqrt((sPri(1)-sPri(2))^2 + (sPri(2)-sPri(3))^2 + (sPri(3)-sPri(1))^2 );

        break;
    end
end

stressOut = [stress(1, 1); stress(2, 2); stress(3, 3); stress(1, 2); stress(2, 3); stress(1, 3)];
end
