%% --- 4. Visualizzazione qualitativa: soluzioni stabili vs instabili ---
%  Non richiesto esplicitamente dal progetto, ma incluso per illustrare
%  visivamente (come Fig. 1/2 del testo) l'effetto della violazione
%  della condizione di stabilita' per theta=0 su mesh diverse.

fprintf('\n');
fprintf('Plot qualitativo: soluzioni stabili vs instabili');
fprintf('\n');

meshCoarse = 'mesh0';   % h grande  -> limite di stabilita' piu' permissivo
meshFine   = 'mesh2';   % h piccolo -> limite di stabilita' piu' severo

% --- calcolo dei dt_stab per le due mesh (come in sezione 3) ---
[~,~,triC,coordC,~,~] = input_data3(meshCoarse);
[BlocC,ClocC,AreaC]   = localBasis3(size(triC,1), triC, coordC);
AC = stiffBuild3(size(triC,1), size(coordC,1), triC, BlocC, ClocC, AreaC);
PC = massBuild3(size(triC,1), size(coordC,1), triC, AreaC);
dtStabCoarse = 2/eigs(AC,PC,1,'lm');

[~,~,triF,coordF,~,~] = input_data3(meshFine);
[BlocF,ClocF,AreaF]   = localBasis3(size(triF,1), triF, coordF);
AF = stiffBuild3(size(triF,1), size(coordF,1), triF, BlocF, ClocF, AreaF);
PF = massBuild3(size(triF,1), size(coordF,1), triF, AreaF);
dtStabFine   = 2/eigs(AF,PF,1,'lm');

fprintf('dt_stab(%s) = %.4e\n', meshCoarse, dtStabCoarse);
fprintf('dt_stab(%s) = %.4e\n', meshFine,   dtStabFine);

% --- caso STABILE: theta=1, dt piccolo (va bene per entrambe le mesh) ---
theta_stable = 1;
dt_stable    = 0.9*min(dtStabCoarse, dtStabFine);   % conservativo

figure('Name','Soluzioni stabili');
[~, ~, uh_c, ~, ~] = solveParabolic(meshCoarse, theta_stable, dt_stable, T, uex, force);
subplot(1,2,1);
trisurf(triC, coordC(:,1), coordC(:,2), uh_c);
shading interp; view(-37.5,30); axis tight;
title(sprintf('%s, \\theta=%.1f (stabile)', meshCoarse, theta_stable));

[~, ~, uh_f, ~, ~] = solveParabolic(meshFine, theta_stable, dt_stable, T, uex, force);
subplot(1,2,2);
trisurf(triF, coordF(:,1), coordF(:,2), uh_f);
shading interp; view(-37.5,30); axis tight;
title(sprintf('%s, \\theta=%.1f (stabile)', meshFine, theta_stable));
sgtitle('Esempio di soluzioni stabili');

% --- caso INSTABILE: theta=0, dt scelto TRA i due limiti di stabilita' ---
%  dt_stab(fine) < dt_unstable < dt_stab(coarse)
%  -> sulla mesh coarse resta (quasi) stabile, sulla mesh fine esplode
theta_unstable = 0;
dt_unstable    = 0.5*(dtStabFine + dtStabCoarse);   % nel mezzo

figure('Name','Soluzioni instabili');
[~, ~, uh_c2, ~, ~] = solveParabolic(meshCoarse, theta_unstable, dt_unstable, T, uex, force);
subplot(1,2,1);
trisurf(triC, coordC(:,1), coordC(:,2), uh_c2);
shading interp; view(-37.5,30); axis tight;
title(sprintf('%s, \\theta=%.1f (dt=%.3e < dt_{stab})', meshCoarse, theta_unstable, dt_unstable));

[~, ~, uh_f2, ~, ~] = solveParabolic(meshFine, theta_unstable, dt_unstable, T, uex, force);
subplot(1,2,2);
trisurf(triF, coordF(:,1), coordF(:,2), uh_f2);
shading interp; view(-37.5,30); axis tight;
title(sprintf('%s, \\theta=%.1f (dt=%.3e > dt_{stab})', meshFine, theta_unstable, dt_unstable));
sgtitle('Esempio di soluzioni instabili');