function [] = plotta3(t,matrice,str)
% Questa funzione realizza il plot di una matrice di tre vettori
% stampandoli in sequenza

% t è il tempo (o l'elemento che finirà sull'asse delle ascisse)

% matrice è la matrice i cui valori finiranno sull'asse delle ordinate

% str è il titolo del grafico

    figure(Name=str)
    subplot(3,1,1);
    plot(t,matrice(:,1),LineWidth=1,Color="r")
    title(str)
    grid
    subplot(3,1,2);
    plot(t,matrice(:,2),LineWidth=1,Color="g")
    grid
    subplot(3,1,3);
    plot(t,matrice(:,3),LineWidth=1,Color="b")
    grid
end
