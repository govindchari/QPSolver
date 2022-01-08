#include "QP.hpp"

QP::QP(Eigen::MatrixXd Qi, Eigen::VectorXd qi, Eigen::MatrixXd Ai, Eigen::VectorXd bi, Eigen::MatrixXd Gi, Eigen::VectorXd hi) : Q(Qi), q(qi), A(Ai), b(bi), G(Gi), h(hi)
    {
        //Variable lengths
        nx = q.size();
        ns = h.size();
        nz = h.size();
        ny = b.size();
        N = nx + ns + nz + ny;

        DELTA delta(nx, ns, nz, ny);

        //KKT System
        KKT.setZero(N, N);
        rhs_a.setZero(N, 1);
        rhs_c.setZero(N, 1);
        p_a.setZero(N, 1);
        p_a.setConstant(2);   //REMOVE THIS LATER      
        p_c.setZero(N, 1);

        //Optimization variables
        x.setZero(nx, 1);
        s.setZero(ns, 1);
        s.setConstant(1);   //REMOVE THIS LATER
        z.setZero(nz, 1);
        y.setZero(ny, 1);

    }

void QP::rhs_kkt_a(){
    Eigen::MatrixXd At = A.transpose();
    Eigen::MatrixXd Gt = G.transpose();
    rhs_a << 
    -(Q*x + q + At*y + Gt*z),
    -(z),
    -(G*x + s - h),
    -(A*x - b);
}
void QP::index_sol_a(){
    delta.x_a = p_a(Eigen::seq(0,nx));
    delta.s_a = p_a(Eigen::seq(nx+1,nx+ns));
    delta.z_a = p_a(Eigen::seq(nx+ns+1,nx+ns+nz));
    delta.y_a = p_a(Eigen::seq(nx+ns+nz+1,N-1));
}
void QP::rhs_kkt_c(double sig, double mu){
    Eigen::VectorXd temp;
    temp.setConstant(ns,1);
    temp = sig*mu*temp;
    rhs_c.setZero(N,1);
    rhs_c(Eigen::seq(nx+1,nx+ns)) = (temp - delta.s_a.cwiseProduct(delta.z_a)).cwiseQuotient(s);
}
void QP::index_sol_c(){
    delta.x_c = p_c(Eigen::seq(0,nx));
    delta.s_c = p_c(Eigen::seq(nx+1,nx+ns));
    delta.z_c = p_c(Eigen::seq(nx+ns+1,nx+ns+nz));
    delta.y_c = p_c(Eigen::seq(nx+ns+nz+1,N-1));
}
double QP::linesearch(Eigen::VectorXd x,Eigen::VectorXd dx){
    return 1.0;
}
void QP::centering_params(double &sig, double &mu){

}
void QP::combine_deltas(){
    delta.x = delta.x_a + delta.x_c;
    delta.s = delta.s_a + delta.s_c;
    delta.z = delta.z_a + delta.z_c;
    delta.y = delta.y_a + delta.y_c;
}
void QP::update_vars(double a){
    x += a * delta.x;
    s += a * delta.s;
    z += a * delta.z;
    y += a * delta.y;
}
void QP::initialize_kkt(){

}
void QP::update_kkt(){

}
void QP::logging(){

}
void QP::solve(){
    index_sol_a();
    rhs_kkt_c(3,2);
    std::cout << rhs_c << std::endl;
}
void QP::solve(bool verbose){

}
