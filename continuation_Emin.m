function continuation_Emin(nelx, nely, volfrac, sd, bc, g, fileID, save_to_folder, same)
    Emins =  [0.5, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0.0001];
    allenergies=[];
    allenergiessame = [];
    %% start continuation
    tic;
    sds = [sd];
    if same
        [y, loop, c, x, energies, energiessame] = topthr(nelx, nely, volfrac, Emins(1), g, sd, bc, 0, 1, fileID, sameEmin);
        allenergiessame(1:length(energiessame)) = energiessame;
    else
        [y, loop, c, x, energies, energiessame] = topthr(nelx, nely, volfrac, Emins(1), g, sd, bc, 0, 1, fileID);
    end
    try
        saveas(gcf, strcat(mydir, meshstr, string(0), '.png'));
    catch err
        disp(err.message);
        save_to_folder = false;
    end
    allenergies(1:loop) = energies;
    numofiterations = loop;
    %%
    for i = 2:length(Emins) 
        if same
            [y, loop, c, x, energies, energiessame]=topthr(nelx, nely, volfrac, Emins(i), g, sd, bc, 1, x, fileID, sameEmin);
            allenergiessame(numofiterations + 1 : numofiterations + length(energiessame)) = energiessame;
        else
            [y, loop, c, x, energies, energiessame]=topthr(nelx, nely, volfrac, Emins(i), g, sd, bc, 1, x, fileID);
        end
        allenergies(numofiterations + 1 : numofiterations + loop) = energies;
        numofiterations = numofiterations + loop;
        if save_to_folder
            saveas(gcf, strcat(mydir, meshstr, string(i-1), '.png'));
        end
        figure;
        plot(numofiterations + 1 : numofiterations + loop, energies, 'o-'); axis tight;
    end  
    toc;
    figure; imshow(1-x); caxis([0,1]); axis equal; axis off; set(gca,'Units','normalized','Position',[0 0 1 1]); 
    if save_to_folder; saveas(gcf, strcat(mydir, meshstr, string(i-1), '_sharp.png')); end
    try
        fprintf(fileID, 'Final Obj: %10.6f', c);
        fclose(fileID);
    catch err
        disp(err.message)
    end
    figure; plot(allenergies, 'o-'); axis tight
    if save_to_folder  
        saveas(gcf, strcat(mydir, meshstr, 'energyplot.png'));
    end
    figure;
    plot(allenergiessame, 'o-'); axis tight
