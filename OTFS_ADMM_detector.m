function x_est =OTFS_ADMM_detector(H,N,M,M_mod,r, rho, alpha,n_ite)
r = reshape(r,N*M,1);
x_est=zeros(M*N,1)+1j*zeros(M*N,1);
x0=zeros(M*N,1)+1j*zeros(M*N,1);
y=zeros(M*N,1)+1j*zeros(M*N,1);   % Lagrange's multiplier
alphabet1 = qammod(0:M_mod-1,M_mod);
alphabet=alphabet1(:);
inv_mat=eye(M*N)/(H'*H+rho*eye(M*N));
Hr=H'*r;
for k=1:n_ite
    term=1/(rho-alpha)*rho*x0+y;
    xq=projection_ADMM(term);
%     x0=inv_mat*(Hr+rho*xq-y);
x0=inv_mat*(Hr+rho*(xq-y));
%     y=y+rho*(x0-xq);
 y=y+rho*(x0+xq);
end
for i=1: length(r)
    x_est(i)=(ml_distance_detection(x0(i),alphabet));
end
end