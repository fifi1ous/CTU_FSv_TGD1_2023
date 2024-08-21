function [sd re] = xy2sd(xy)
% --------------------------------------------------------------------
% prevede rovinne kartezske souradnice [X Y] pomoci rovnic Krovakova
% zobrazeni na kartograficke souradnice [S D] na kouli
% postup transformaci a zobrazeni: ([X Y]->[ro epsilon]->[S D])
% 
% IN:  xy = [X Y] ... matice rovinnych kartezskych souradnic X[m], Y[m]
%                    (na radcich postupne X Y)
%                     
% OUT: sd = [S D]        ... matice kartografickych souradnic
%                            kart.sirka S[rad], kart.delka D[rad]
%                           (na radcich postupne S D)
%      re = [ro epsilon] ... matice rovinnych polarnich souradnic
%                            pruvodic ro[m], polarni uhel epsilon[rad]
%                           (na radcich postupne ro epsilon)
% --------------------------------------------------------------------
% functions used:
% rad2dms
% --------------------------------------------------------------------

% kontrola poctu parametru
if (nargin ~= 1)
   error('xy2sd: Chybny pocet vstupnich parametru');
end
% kontrola rozmeru
[k,m]=size(xy);
if (m ~= 2)
   error('xy2sd: Vstupni matice musi obsahovat souradnice X a Y');
end

% definice konstant
S0  = 78.5*pi/180;
Rk = 6380703.6105;
Rc = 0.9999 * Rk;
ro0 = Rc / tan(S0);

% prevod na polarni souradnice (ro, epsilon)
ro      = sqrt( xy(:,1).^2 + xy(:,2).^2 );
epsilon = atan( xy(:,2) ./ xy(:,1) );
          % pozn.: vzdy muze nastat jen 1.kvadrant

% prevod na kartograficke souradnice (S, D)
D = epsilon ./ sin(S0);
argument = (ro0./ro).^(1/sin(S0)) .* tan(S0/2+pi/4);
S = 2*( atan(argument) - pi/4 );

% vystupni matice
sd = [S D];
sd=sd/pi*180;
re = [ro epsilon/pi*180];
end
