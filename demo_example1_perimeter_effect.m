% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID=-1;
% whether calculate compliance with respect to the same Emin
sameEmin = 5e-5;
same = false;
% parameters
volfrac = 0.5;
gs = [0, 1e-7, 1.4e-7, 6e-7];
nelx = 400;
nely = 200;
bc = 'cantilever_rb';
sd = 1;
for g = gs
    continuation_Emin(nelx, nely, volfrac, sd, bc, g, fileID, save_to_folder, same)
end

