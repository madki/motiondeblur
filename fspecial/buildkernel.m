% function to build kernel
% K = camera intrinsic matrix 3x3
% R rotation matrix 3x3
% T translation vector 3x1
% d scene depth, scalar
% N normal to image plane 3x1
% nr = number of rows in image
% nc = number of columns in image
function [A] = buildkernel(K, R, T, d, N, nr, nc)
    A = zeros(nr*nc);
    H = K*(R+T*N'/d)/K;
    for j=1:nc
        for i=1:nr
            temp = H\[i, j, 1]';
            tempi = int64(temp);
            if temp(1)>=tempi(1)
                ul = tempi(1);
                uh = ul + 1;
            else
                uh = tempi(1);
                ul = uh -1;
            end
            if temp(2)>=tempi(2)
                vl = tempi(2);
                vh = vl + 1;
            else
                vh = tempi(2);
                vl = vh -1;
            end
            A(i*j, ul*vl) = (uh - temp(1))*(vh - temp(1));
            A(i*j, uh*vl) = (temp(1) - ul)*(vh - temp(1));
            A(i*j, uh*vh) = (temp(1) - ul)*(temp(1) - vl);
            A(i*j, ul*vh) = (uh - temp(1))*(temp(1) - vl);
        end
    end
end