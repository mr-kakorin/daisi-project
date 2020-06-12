R = table2array(readtable("R.txt"));
Phi = table2array(readtable("Phi.txt"));
T = table2array(readtable("T.txt"));
R = reshape(R, 240, 240);
Phi = reshape(Phi, 240, 240);
T = reshape(T, 240, 240);
pcolor(R,Phi,T);
shading interp
colorbar