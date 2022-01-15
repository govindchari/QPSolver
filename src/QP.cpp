#include "QP.hpp"
#include <iostream>
#include <Eigen/SparseCholesky>

QP::QP(Eigen::MatrixXd Qi, Eigen::VectorXd qi, Eigen::MatrixXd Ai, Eigen::VectorXd bi, Eigen::MatrixXd Gi, Eigen::VectorXd hi) : Q(Qi), q(qi), A(Ai), b(bi), G(Gi), h(hi)
{

    //Variable lengths
    nx = q.size();
    ns = h.size();
    nz = h.size();
    ny = b.size();
    N = nx + ns + nz + ny;

    DELTA delta(nx, ns, nz, ny);

    //Regularization Matrix
    double eps = 1e-7;
    Eigen::MatrixXd reg_dense(N, N);
    reg_dense.setZero(N, N);
    reg_dense(Eigen::seq(0, nx + ns - 1), Eigen::seq(0, nx + ns - 1)) = eps * Eigen::MatrixXd::Identity(nx + ns, nx + ns);
    reg_dense(Eigen::seq(nx + ns, N - 1), Eigen::seq(nx + ns, N - 1)) = -eps * Eigen::MatrixXd::Identity(nz + ny, nz + ny);
    reg = reg_dense.sparseView();

    //KKT System
    KKT.setZero(N, N);
    rhs_a.setZero(N, 1);
    rhs_c.setZero(N, 1);
    p_a.setZero(N, 1);
    p_c.setZero(N, 1);

    //Optimization variables
    x.setZero(nx, 1);
    s.setZero(ns, 1);
    z.setZero(nz, 1);
    y.setZero(ny, 1);
}

void QP::rhs_kkt_a()
{
    // Unconstrained Optimization
    if (A.size() == 0 && G.size() == 0)
    {
        rhs_a << -(Q.sparseView() * x + q);
    }
    // Inequality Constrained Optimization
    else if (A.size() == 0 && G.size() != 0)
    {
        rhs_a << -(Q.sparseView() * x + q + (G.transpose()).sparseView() * z),
            -(z),
            -(G.sparseView() * x + s - h);
    }
    // Equality Constrained Optimization
    else if (A.size() != 0 && G.size() == 0)
    {
        rhs_a << -(Q * x + q + (A.transpose()).sparseView() * y),
            -(A.sparseView() * x - b);
    }
    // Equality and Inequality Constrained Optimization
    else
    {
        rhs_a << -(Q * x + q + A.transpose().sparseView() * y + G.transpose().sparseView() * z),
            -(z),
            -(G.sparseView() * x + s - h),
            -(A.sparseView() * x - b);
    }
}
void QP::index_sol_a()
{
    delta.x_a = p_a(Eigen::seq(0, nx - 1));
    delta.s_a = p_a(Eigen::seq(nx, nx + ns - 1));
    delta.z_a = p_a(Eigen::seq(nx + ns, nx + ns + nz - 1));
    delta.y_a = p_a(Eigen::seq(nx + ns + nz, N - 1));
}
void QP::rhs_kkt_c(const double sig, const double mu)
{
    Eigen::VectorXd temp;
    temp.setConstant(ns, 1);
    temp = sig * mu * temp;
    rhs_c.setZero(N, 1);
    if (G.size() != 0)
    {
        rhs_c(Eigen::seq(nx + 1, nx + ns)) = (temp - delta.s_a.cwiseProduct(delta.z_a)).cwiseQuotient(s);
    }
}
void QP::index_sol_c()
{
    delta.x_c = p_c(Eigen::seq(0, nx - 1));
    delta.s_c = p_c(Eigen::seq(nx, nx + ns - 1));
    delta.z_c = p_c(Eigen::seq(nx + ns, nx + ns + nz - 1));
    delta.y_c = p_c(Eigen::seq(nx + ns + nz, N - 1));
}
double QP::linesearch(const Eigen::VectorXd x, const Eigen::VectorXd dx)
{
    int len = x.size();
    Eigen::VectorXd temp(len);
    for (int i = 0; i < len; i++)
    {
        if (dx[i] < 0)
        {
            temp[i] = -x[i] / dx[i];
        }
        else
        {
            temp[i] = std::numeric_limits<double>::max();
        }
    }
    return fmin(1.0, temp.minCoeff());
}
void QP::centering_params(double &sig, double &mu)
{
    mu = s.dot(z) / ns;
    double a = fmin(linesearch(s, delta.s_a), linesearch(z, delta.z_a));
    sig = std::pow((s + a * delta.s_a).dot(z + a * delta.z_a) / s.dot(z), 3);
}
void QP::combine_deltas()
{
    delta.x = delta.x_a + delta.x_c;
    delta.s = delta.s_a + delta.s_c;
    delta.z = delta.z_a + delta.z_c;
    delta.y = delta.y_a + delta.y_c;
}
void QP::update_vars(double a)
{
    x += a * delta.x;
    s += a * delta.s;
    z += a * delta.z;
    y += a * delta.y;
}
void QP::initialize_kkt()
{
    auto idx_x = Eigen::seq(0, nx - 1);
    auto idx_s = Eigen::seq(nx, nx + ns - 1);
    auto idx_z = Eigen::seq(nx + ns, nx + ns + nz - 1);
    auto idx_y = Eigen::seq(nx + ns + nz, N - 1);

    KKT(idx_x, idx_x) = Q;

    if (G.size() != 0)
    {
        KKT(idx_x, idx_z) = G.transpose();
        KKT(idx_z, idx_x) = G;
        KKT(idx_z, idx_s) = Eigen::MatrixXd::Identity(nz, nz);
        KKT(idx_s, idx_z) = Eigen::MatrixXd::Identity(ns, ns);
    }
    if (A.size() != 0)
    {
        KKT(idx_x, idx_y) = A.transpose();
        KKT(idx_y, idx_x) = A;
    }
}
void QP::update_kkt()
{
    if (G.size() != 0)
    {
        auto idx_s = Eigen::seq(nx, nx + ns - 1);
        Eigen::VectorXd temp = z.cwiseQuotient(s);
        KKT(idx_s, idx_s) = temp.asDiagonal();
    }
}
void QP::regularize_kkt()
{
    KKT_reg = KKT.sparseView();
    KKT_reg += reg;
}
void QP::iterative_refinement(Eigen::VectorXd &sol, const Eigen::VectorXd &rhs, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver)
{
    Eigen::VectorXd dl(rhs.size());
    while ((KKT * sol - rhs).norm() > tol.iter_ref)
    {
        dl = solver.solve((rhs - KKT * sol).sparseView());
        sol += dl;
    }
}
void QP::logging(const int iter, const double a)
{
    double temp1 = x.transpose() * Q * x;
    double temp2 = q.transpose() * x;
    auto J = temp1 + temp2;
    double eq_res = 0;
    double ineq_res = 0;
    double gap = 0;

    if (A.size() != 0)
    {
        eq_res = (A * x - b).norm();
    }
    if (G.size() != 0)
    {
        gap = s.transpose() * z;
        ineq_res = (G * x + s - h).norm();
    }

    printf("%3d    %10.3e    %9.2e    %9.2e    %9.2e    %6.4f\n", iter, J, gap, eq_res, ineq_res, a);
}
void QP::initialize()
{
    auto idx_x = Eigen::seq(0, nx - 1);
    auto idx_s = Eigen::seq(nx, nx + ns - 1);
    auto idx_y = Eigen::seq(nx + ns, nx + ns + ny - 1);
    auto Ni = nx + nz + ny;

    Eigen::MatrixXd A_init;
    Eigen::VectorXd rhs(Ni);

    A_init.setZero(Ni, Ni);
    A_init(idx_x, idx_x) = Q;

    if (G.size() != 0)
    {
        A_init(idx_x, idx_s) = G.transpose();
        A_init(idx_s, idx_x) = G;
        A_init(idx_s, idx_s) = -Eigen::MatrixXd::Identity(ns, ns);
    }

    if (A.size() != 0)
    {
        A_init(idx_x, idx_y) = A.transpose();
        A_init(idx_y, idx_x) = A;
    }

    rhs << -q, h, b;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_init.sparseView());
    Eigen::VectorXd sol = solver.solve(rhs);

    x = sol(idx_x);
    z = sol(idx_s);
    y = sol(idx_y);

    auto neg_z = -z;
    auto a_p = -(neg_z.minCoeff());

    Eigen::VectorXd ones(ns);
    ones.setConstant(1);

    if (a_p < 0)
    {
        s = -z;
    }
    else
    {
        s = -z + (1 + a_p) * ones;
    }

    auto a_d = -z.minCoeff();
    if (a_d >= 0)
    {
        z += (1 + a_d) * ones;
    }
}
void QP::solve()
{
    solve(true);
}
void QP::solve(bool verbose)
{
    if (verbose)
    {
        printf("iter      objv          gap         |Ax-b|      |Gx+s-h|     step\n");
        printf("------------------------------------------------------------------\n");
    }

    //Stopping criterion stuff
    unsigned int iter = 1;
    double Jprev = std::numeric_limits<double>::max();
    double temp1;
    double temp2;
    double Jcurr;
    double ineq_res = 0;
    double eq_res = 0;
    double gap = 0;
    bool proceed = false;

    initialize();
    initialize_kkt();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    double sig;
    double mu;
    double a;
    do
    {
        update_kkt();
        regularize_kkt();

        //Affine scaling step
        rhs_kkt_a();
        solver.compute(KKT_reg);
        p_a = solver.solve(rhs_a);
        iterative_refinement(p_a, rhs_a, solver);
        index_sol_a();

        //Centering and correction step
        centering_params(sig, mu);
        rhs_kkt_c(sig, mu);
        p_c = solver.solve(rhs_c);
        iterative_refinement(p_c, rhs_c, solver);
        index_sol_c();

        //Combine prediction and correction steps
        combine_deltas();

        //Update solution iterate after linesearch
        a = fmin(1, 0.99 * fmin(linesearch(s, delta.s), linesearch(z, delta.z)));
        update_vars(a);

        if (verbose)
        {
            logging(iter, a);
        }

        //Stopping criterion
        temp1 = x.transpose() * Q * x;
        temp2 = q.transpose() * x;
        Jcurr = temp1 + temp2;
        if (A.size() != 0)
        {
            eq_res = (A * x - b).norm();
        }
        if (G.size() != 0)
        {
            gap = s.transpose() * z;
            ineq_res = (G * x + s - h).norm();
        }

        //Evaluate stopping criteria
        proceed = ((std::abs(Jcurr - Jprev) > tol.cost) || (eq_res > tol.constraint) || (ineq_res > tol.constraint) || (gap > tol.gap)) && (iter < 25);

        Jprev = Jcurr;
        iter++;
    } while (proceed);

    //Store solution
    solution.x = x;
    solution.s = s;
    solution.z = z;
    solution.y = y;
}
