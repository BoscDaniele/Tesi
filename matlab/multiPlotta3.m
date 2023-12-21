function [] = multiPlotta3(t,matrice1, matrice2, str1, str2)
% Questa funzione realizza il plot di due matrici di tre vettori
% stampandoli in sequenza una accanto all'altra

% t è il tempo (o l'elemento che finirà sull'asse delle ascisse)

% matrice 1 e 2 sono le matrici i cui valori finiranno sugli assi delle
% ordinate

% str 1 e 2 sono i titoli dei grafici


    figure(Name=str1)
    subplot(3,2,1);
    plot(t,matrice1(:,1),LineWidth=1,Color="r")
    title(str1)
    grid
    subplot(3,2,3);
    plot(t,matrice1(:,2),LineWidth=1,Color="g")
    grid
    subplot(3,2,5);
    plot(t,matrice1(:,3),LineWidth=1,Color="b")
    grid
    subplot(3,2,2);
    plot(t,matrice2(:,1),LineWidth=1,Color="r")
    title(str2)
    grid
    subplot(3,2,4);
    plot(t,matrice2(:,2),LineWidth=1,Color="g")
    grid
    subplot(3,2,6);
    plot(t,matrice2(:,3),LineWidth=1,Color="b")
    grid
end