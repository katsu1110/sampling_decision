function E = discretize_signal(E)

nX =  size(E.X, 2);
% O_pref=1;
ixp=1; ixa=1+nX/2;
pall = squeeze(E.Signal(:,ixp,:));
pall = pall(:);
nall = squeeze(E.Signal(:,ixa,:));
nall = nall(:);
cutval = 0.001*[-4, -1, 0, 1.2, 1.9, 2.7, 3.4, 4.7, 5.3, 8];
binval = [-pi/2, -3*pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2, 7*pi/8, pi];
pall = arrayfun(@(x) binnize(x, cutval, binval), pall);
nall = arrayfun(@(x) binnize(x, cutval, binval), nall);
E.Signal(:,ixp,:) = reshape(pall, size(E.Signal(:,ixp,:)));
E.Signal(:,ixa,:) = reshape(nall, size(E.Signal(:,ixa,:)));


function a = binnize(a,cutval, binval)
if a <= cutval(2)
    a = binval(1);
elseif a > cutval(2) && a <= cutval(3)
    a = binval(2);
elseif a > cutval(3) && a <= cutval(4)
    a = binval(3);
elseif a > cutval(4) && a <= cutval(5)
    a = binval(4);
elseif a > cutval(5) && a <= cutval(6)
    a = binval(5);
elseif a > cutval(6) && a <= cutval(7)
    a = binval(6);
elseif a > cutval(7) && a <= cutval(8)
    a = binval(7);
elseif a > cutval(8) && a <= cutval(9)
    a = binval(8);
    elseif a > cutval(9)
    a = binval(9);
end
