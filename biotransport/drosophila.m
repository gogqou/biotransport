function res = drosophila()
    clear
    clf
    D = 10;
    a = 1;
%     tau = 5000;
    N = 20;
    M =25;
    ubound = 5;
    umin = .1;
    vbound = 2*pi;
    vmin = .1;
    du = (ubound-umin)/(N+1);
    dv = (vbound-vmin)/(M+1);
    u = umin:du:ubound;
    v = vmin:dv:vbound;
    Tstart = 0;           % Start time (time) 
    Tend=10;              % End time (time)
    dt = 1;
    C0 = zeros(N,M);  % Define the initial condition column vector
%     C0(N,M) = 1;          % Initial conditions 
    C0(N,round(M/2)) = 1;
    init = reshape(C0, M*N, 1); %reshape into a column vector
%     init(N*M)=1;
    lapconst= zeros(N,M); %allocate laplacian constant (the 1/ a^2 sinhu ....part)
    x_matrix = zeros(N,M); %allocate 
    y_matrix= zeros(N,M); %allocate
    for i = 1:N
        for j = 1:M
            vsq=sin(v(j))^2;
            usq=sinh(u(i))^2;
            lapconst(i,j) = 1/(a^2*(vsq+usq));
            x_matrix(i,j) = a*cosh(u(i))*cos(v(j)); %convert to cartesian
            y_matrix(i,j) = a*sinh(u(i))*sin(v(j)); %convert to cartesian
        end
    end
%     lapconst
%     x_vector = reshape(x_matrix, M*N,1);
%     y_vector = reshape(y_matrix, M*N,1);
%     x = [-x_matrix x_matrix];
%     y = [-y_matrix y_matrix];
%     [x_vector y_vector]
    [t,Y] = ode45(@odes,[Tstart:dt:Tend],init);

    for i = 1:length(t)
        C = reshape(Y(i,:),N,M);
        surf(x_matrix,y_matrix,C);
    end

    set(0, 'defaultaxesfontsize', 14);
    xlabel ('X', 'FontSize', 15);
    ylabel ('Y', 'FontSize', 15);
    title('Dextran', 'FontSize', 16);
%     initial = reshape(Y(1,:), N,M)
    res = reshape(Y(end,:), N,M);
    
function derivs = odes(t,Y)
%     t
    uv = reshape(Y,N,M);
    uv(N,round(M/2)) = 1;
    duvdt  = zeros(N,M);
    duvdu2 = zeros(N,M);
    duvdv2 = zeros(N,M);
    
    duvdu2(1,:) = (uv(2,:) - uv(1,:))/(du^2); 
    duvdu2(N,:) = (-uv(N,:) +uv(N-1,:))/(du^2);
    
    duvdv2(:,1) = (uv(:,2) - uv(:,1))/(dv^2);
    duvdv2(:,end) = (-uv(:,end)+ uv(:,end-1))/(dv^2);

    duvdu2 (2:N-1,:) = (uv(3:N,:)-2*uv(2:N-1,:) + uv(1:N-2,:))/(du^2);
    duvdv2 (:,2:M-1) = (uv(:,3:M)-2*uv(:,2:M-1) + uv(:,1:M-2))/(dv^2);
    
%     reactnterm = 1/tau* uv;
%     duvdt = D*lapconst.*(duvdu2 + duvdv2)-reactnterm;
    duvdt = D*lapconst.*(duvdu2 + duvdv2);
%     duvdt (N, round(M/2)) = 1;
%     duvdt
    derivs = reshape(duvdt, M*N, 1); 
end

end