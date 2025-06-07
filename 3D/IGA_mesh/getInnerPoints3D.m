function innerPoints = getInnerPoints3D(geomData, numInner)
% Generate inner points for paper1.
%

shapeType = geomData{1};
shapeParams = geomData{2};

switch shapeType
    case 'cub'
        a = shapeParams(4);
        b = shapeParams(5);
        c = shapeParams(6);
        O = shapeParams(1 : 3);
        nX = numInner(1);
        nY = numInner(2);
        nZ = numInner(3);
        innerPoints = zeros(nX*nY*nZ, 3);
        count = 1;
        for ix = 1 : nX
            for iy = 1 : nY
                for iz = 1 : nZ
                    x = (ix - 0.5) * 2 * a / nX - a;
                    y = (iy - 0.5) * 2 * b / nY - b;
                    z = (iz - 0.5) * 2 * c / nZ - c;
                    innerPoints(count, :) = [x, y, z];
                    count = count + 1;
                end
            end
        end
        innerPoints = innerPoints + O;
    case 'star_e3'
        O = shapeParams(1 : 3);
        [H1, H2, H3, R1, R2] = deal(shapeParams(4), shapeParams(5), ...
            shapeParams(6), shapeParams(7), shapeParams(8));
        nR1 = numInner(1);
        nR2 = numInner(2);
        nT = numInner(3);
        nY = numInner(4);
        nZ = numInner(5);
        RMid = R1 - H2;

        innerPoints = [];

        for k = 1 : nZ
            z = k * H3 / (nZ + 1);
            for i = 1 : nR2
                for j = 1 : nY
                    y = j * R2 / nY;
                    x = R1 - (sqrt(R1^2 - y^2) - sqrt(R2^2 - y^2) + R2 - RMid) / (nR2 + 1) * i;
                    x1 = x;
                    y1 = y;
                    [innerPoints] = [innerPoints; x y z; y1 x1 z];
                end
            end

            t1 = 0;
            t2 = pi/2-t1;

            for i = 1 : nR2 + nR1
                for j = 1 : nT - 1
                    t = j * (t2 - t1) / nT;
                    d1 = sqrt(2) * R2 * cos(pi/4 - t) - R2 * sqrt(cos(pi/2 - 2 * t));
                    d2 = sqrt(2) * R2 * cos(pi*3/4 + t) + sqrt(R2^2 * cos(pi*3/2 + 2*t) + R1^2);
                    x2 = ((d2 - d1) * i / (nR2 + nR1 + 1) + d1) * cos(t) + R2;
                    y2 = ((d2 - d1) * i / (nR2 + nR1 + 1) + d1) * sin(t) + R2;
                    [innerPoints] = [innerPoints; x2 y2 z];
                end
            end
        end
        innerPoints = innerPoints + O;
    case 'star_e4'
        O = shapeParams(1 : 3);
        [H1, H2, H3, R1, R2] = deal(shapeParams(4), shapeParams(5), ...
            shapeParams(6), shapeParams(7), shapeParams(8));
        nR1 = numInner(1);
        nR2 = numInner(2);
        nT = numInner(3);
        nY = numInner(4);
        nZ = numInner(5);
        RMid = R1 - H2;

        alph = 72;

        innerPoints = [];

        for k = 1 : nZ
            z = k * H3 / (nZ + 1);
            for i = 1 : nR2
                for j = 1 : nY
                    y = j * R2 / nY;
                    x = R1 - (sqrt(R1^2 - y^2) - sqrt(R2^2 - y^2) + R2 - RMid) / (nR2 + 1) * i;
                    x1 = y * cosd(90 - alph) + x * sind(90 - alph);
                    y1 = x * cosd(90 - alph) - y * sind(90 - alph);
                    [innerPoints] = [innerPoints; x y z; x1 y1 z];
                end
            end

            t1 = 0;
            t2 = alph-t1;

            for i = 1 : nR2 + nR1
                for j = 1 : nT - 1
                    t = j * (t2 - t1) / nT;
                    d_1 = R2 / sind(alph / 2);
                    d1 = d_1 * cosd(alph / 2 - t) - sqrt(d_1^2 * ((cosd(alph/2 - t))^2 - 1) + R2^2);
                    d2 = d_1 * cosd(180 - alph / 2 + t) + sqrt(d_1^2 * ((cosd(180 - alph/2 + t))^2 - 1) + R1^2);
                    x2 = ((d2 - d1) * i / (nR2 + nR1 + 1) + d1) * cosd(t) + R2 / tand(alph / 2);
                    y2 = ((d2 - d1) * i / (nR2 + nR1 + 1) + d1) * sind(t) + R2;
                    [innerPoints] = [innerPoints; x2 y2 z];
                end
            end
        end
        innerPoints = innerPoints + O;
end

end

