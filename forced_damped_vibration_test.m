m = 12.5;
c = 0.01;
k = 4;
F0 = 3.75;
omega_0 = 6.2;

omega_n = sqrt(k/m);
c_cr = 2*m*omega_n;

A = F0/(-m*omega_0^2-((c*omega_0)^2/(m*omega_0^2-k))+k);
B = ((c*omega_0)/(m*omega_0^2-k))*A;

D1 = sqrt(A^2+B^2)
phi1 = atan(B/A)

D2 = (F0/k)/sqrt((1-(omega_0/omega_n)^2)^2+(2*(c/c_cr)*(omega_0/omega_n))^2)
phi2 = atan( (2*(c/c_cr)*(omega_0/omega_n)) / (1-(omega_0/omega_n)^2) )