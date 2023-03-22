% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID=-1;
% whether calculate compliance with respect to the same Emin
sameEmin = 5e-5;
same = false;
% parameters
volfrac = 0.5;
nelx = 200;
nely = 100;
bc = 'cantilever_rb';
g = 0;
sd = 1;
continuation_Emin(nelx, nely, volfrac, sd, bc, g, fileID, save_to_folder, same)