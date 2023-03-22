%%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
function [y, loop, c,  x, energies, energiessame]=topthr(nelx, nely, volfrac, frac, g, sd, bc, continuation, x, fileID, sameEmin)
% nelx: number of elements on x axis
% nely: number of elements on y axis
% volfrac: volume fraction of material to total area
% frac: fraction of Emin to E0
% sd1: filter size for thresholding
% sd2: filter size for solving equation
% sd3: filter size for the perimeter term
% bc:  takes values in {'cantilever_rb', 'cantilever_central', 'mbb'},
%      enforces BCs
% g: coefficient of the perimeter term
% continuation: 1 indicates to use x as the initial guess
%               0 indicates to use constant density as the initial guess
% x: initial guess, if continuaton is 0, can be anything.
% fileID: the opened file to log outputs
% If you don't want to save the results, just let fileID to be negative 
% and the results will be displayed in the command window
% sameEmin: calculate the compliance each iteration using this specified value as Emin
% return c: compliance computed with filtered chi
%% MATERIAL PROPERTIES
E0 = 5000*8/3; 
Emin = E0*frac;
nu = 1/3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = (1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]));
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
switch bc
    case 'mbb'
        F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1); 
        fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
    case 'cantilever_rb'
        force_length = floor(nelx/8); F = sparse(2*(nely+1)*((nelx+1-force_length):(nelx+1)), ones(force_length+1, 1), -250/nely*ones(force_length+1, 1), 2*(nely+1)*(nelx+1), 1);% 
        fixeddofs = [1 : 2*(nely+1)]; 
    case 'cantilever_central'
        F = sparse(2*(nely + 1) * nelx + 2*(floor(nely/2) + 1), 1, -1, 2*(nely + 1)*(nelx + 1),1);
        fixeddofs = [1 : 2*(nely+1)];  
end
U = zeros(2*(nely + 1)*(nelx + 1), 1); 
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs); 
%% INITIALIZE ITERATION
if continuation == 0
    x = repmat(volfrac,nely,nelx);
    xPhys = x;
else
    if sd > 0
        xPhys = imgaussfilt(x, sd);
    else
        xPhys = x;
    end
end
M = floor(nelx*nely*volfrac);
loop = 0;
change = 100;
gamma = g*sqrt(2*pi)/(sd*1/nely);
energies = [];
energiessame = [];
%% START ITERATION
figure('Renderer', 'painters', 'Position', [10 10 100 nely/nelx*100]);
print_to_file= true;
try 
    fprintf(fileID, 'Displaying\n');
catch err
    disp(err.message)
    disp('Now display the output in the command window.')
    print_to_file=false;
end
while change>0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*1./(1/Emin+xPhys(:)'*(1/E0-1/Emin)), 64*nelx*nely, 1);
  K = sparse(iK,jK,sK); K = (K+K')/2; 
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);     
  c = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ce));
  energies(loop) = c;
  %% threshold dynamics--calculate c sigma:sigma
  csigmasigmae = 1./(1/Emin+xPhys*(1/E0-1/Emin)).^2.*ce;
  %% threshold dynamics--filtering
  if sd > 0
    csigmasigmae = imgaussfilt(csigmasigmae, sd, 'Padding', 'symmetric');
  end
  perterm = 1-2*x;
  if sd > 0 
    perterm = imgaussfilt(perterm, sd);
  end
  %% threshold dynamics--thresholding
  phi = (1/E0 - 1/Emin)*csigmasigmae + gamma*perterm;
  [sortedphi,I] = sort(phi(:));
  thr = sortedphi(M);
  % thresholding
  if thr > 0   
    xnew = 1 - double(phi>0);
  else
    xnew = zeros(nelx*nely, 1);
    xnew(I(1:M)) = 1;
    xnew = reshape(xnew, nely, nelx);
  end
  if nargin > 10
    sK = reshape(KE(:)*((sameEmin+xnew(:)'*(E0 - sameEmin))), 64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2; 
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); 
    energiessame(loop) = sum(sum((sameEmin+xnew*(E0-sameEmin)).*ce));   
  end
  if sd > 0
    xPhys = imgaussfilt(xnew, sd, 'Padding', 'symmetric');
  else
    xPhys = xnew;
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  if print_to_file
    fprintf(fileID, ' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change);
  else
    fprintf(' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change);
  end
  %% PLOT DENSITIES
  colormap(gray); imshow(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
sK = reshape(KE(:)*1./(1/Emin+xPhys(:)'*(1/E0-1/Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2; 
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);     
c = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ce));
loop = loop + 1;
energies(loop) = c;
if nargin > 10
    sK = reshape(KE(:)*((sameEmin+xnew(:)'*(E0 - sameEmin))), 64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2; 
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); 
    energiessame(loop) = sum(sum((sameEmin+xnew*(E0-sameEmin)).*ce));   
end
sK = reshape(KE(:)*(x(:)'*(E0-Emin) + Emin), 64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2; 
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%% FINAL OBJECTIVE FUNCTION WITHOUT SMOOTHING
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); 
y = sum(sum(((E0-Emin)*x+Emin).*ce));
if print_to_file
    fprintf(fileID, ' sharp interface energy: %11.6f\n', y);
else
    fprintf(' sharp interface energy: %11.8f\n', y);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written based on 88 line MATLAB code written by 
% E. Andreassen, A. Clausen, M. Schevenels,
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
% Technical University of Denmark,                                         %
% DK-2800 Lyngby, Denmark.                                                 %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%