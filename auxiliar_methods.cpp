
// ----
// 2024 - 01 - 10 --> Para algoritmo POT version actualizada
#include "pot-library-cpp-wrapped/network_simplex_simple.h"
#include "pot-library-cpp-wrapped/network_simplex_simple_omp.h"
#include "pot-library-cpp-wrapped/EMD.h"
#include <cstdint>
// ----------------------

// --------------------------------------------------
// Methods for Optimal Transport
// 1. Method for the cost
double cost( std::vector<double> M,  std::vector<double> G, int n, double reg){
    // multiply M and G
    // sum all elements of the multiplication of M and G
    double sum = 0;
    for(int i=0; i<n*n; i++){
        sum += M[i] * G[i];
    }

    // apply (1/2)x² (element-wise) to all elements of G and sum them
    double sum2 = 0;
    for(int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            sum2 += G[i*n+j]*G[i*n+j];
        }
    }
    sum2 = sum2/2;

    // return the sum of the two sums
    return sum + reg*sum2;
}
// 2. Method EMD 
int EMD( int n1, int n2, std::vector<double> X, std::vector<double> Y, std::vector<double> D, std::vector<double> &G, std::vector<double> alpha, std::vector<double> beta, double cost, uint64_t maxIter)  {

    using namespace lemon;
    uint64_t n, m, cur;
    typedef FullBipartiteDigraph Digraph;
    DIGRAPH_TYPEDEFS(Digraph);
    // Get the number of non zero coordinates for r and c
    n=0;
    for (int i=0; i<n1; i++) {
        //double val=*(X+i);
        double val=X[i];
        if (val>0) {   
            n++;
        }else if(val<0){
			return INFEASIBLE;
		}
    }
    m=0;
    for (int i=0; i<n2; i++) {
        //double val=*(Y+i);
        double val=Y[i];
        if (val>0) {
            m++;
        }else if(val<0){
			return INFEASIBLE;
		}
    }
    // Define the graph
    std::vector<uint64_t> indI(n), indJ(m);
    std::vector<double> weights1(n), weights2(m);
    Digraph di(n, m);
    NetworkSimplexSimple<Digraph,double,double, node_id_type> net(di, true, (int) (n + m), n * m, maxIter);
    // Set supply and demand, don't account for 0 values (faster)
    cur=0;
    for (uint64_t i=0; i< (uint64_t) n1; i++) {
        double val=X[i];
        if (val>0) {
            weights1[ cur ] = val;
            indI[cur++]=i;
        }
    }
    // Demand is actually negative supply...
    cur=0;
    for (uint64_t i=0; i< (uint64_t) n2; i++) {
        double val=Y[i];
        if (val>0) {
            weights2[ cur ] = -val;
            indJ[cur++]=i;
        }
    }
    net.supplyMap(&weights1[0], (int) n, &weights2[0], (int) m);
    // Set the cost of each edge
    int64_t idarc = 0;
    for (uint64_t i=0; i<n; i++) {
        for (uint64_t j=0; j<m; j++) {
            double val=D[indI[i]*n2+indJ[j]];
            net.setCost(di.arcFromId(idarc), val);
            ++idarc;
        }
    }
    // Solve the problem with the network simplex algorithm
    int ret=net.run();
    uint64_t i, j;
    if (ret==(int)net.OPTIMAL || ret==(int)net.MAX_ITER_REACHED) {
        //*cost = 0;
        cost = 0;
        Arc a; di.first(a);
        for (; a != INVALID; di.next(a)) {
            i = di.source(a);
            j = di.target(a);
            double flow = net.flow(a);
            cost += flow * (D[indI[i]*n2+indJ[j-n]]);
            G[indI[i]*n2+indJ[j-n]] = flow;
            alpha[indI[i]] = -net.potential(i);
            beta[indJ[j-n]] = net.potential(j);
        }
    }
    return ret;
}

// 3. Method emd
std::vector<double> emd( std::vector<double> a,  std::vector<double> b, std::vector<double> M){
    int n = a.size();
    std::vector<double> alpha(n);
    std::vector<double> beta(n);
    std::vector<double> G((int)(n*n));
    double cost = 0;
    int ret = EMD(n, n, a, b, M, G, alpha, beta, cost, 100000);
    return G;
}
// 4. Method phi 
double phi( int n,  std::vector<double> xk,  std::vector<double> pk, double alpha){
    double sum = 0;
    if (alpha < 0){
        for (int i=0; i < n*n; i++){
            double val = xk[i];
            sum += val * val;
        }
    }
    else {
        for (int i=0; i < n*n; i++){
            double val = xk[i] + alpha * pk[i];
            sum += val * val;
        }
    }
    return sum/2;
}
// 5. Method scalar_search_armijo
std::pair<double, double> scalar_search_armijo( int n, std::vector<double> xk, std::vector<double> pk, double phi0, double derphi0){
    double c1 = 1e-4;
    double alpha0 = 0.99;
    double amin = 0;

    double phi_a0 = phi(n, xk, pk, alpha0);

    // check the break condition
    if (phi_a0 <= phi0 + c1*alpha0*derphi0){
        return std::make_pair(alpha0, phi_a0);
    }
    // otherwise, compute the minimizer of the quadratic interpolant
    double alpha1 = -(derphi0) * (alpha0*alpha0) / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0);
    double phi_a1 = phi(n, xk, pk, alpha1);

    // check the break condition
    if (phi_a1 <= phi0 + c1*alpha1*derphi0){
        return std::make_pair(alpha1, phi_a1);
    }

    while (alpha1 > amin){

        double factor = alpha0*alpha0*alpha1*alpha1*(alpha1-alpha0);
        double a = alpha0*alpha0*(phi_a1-phi0-derphi0*alpha1) - alpha1*alpha1*(phi_a0-phi0-derphi0*alpha0);
        a = a / factor;
        double b = -alpha0*alpha0*alpha0*(phi_a1-phi0-derphi0*alpha1) + alpha1*alpha1*alpha1*(phi_a0-phi0-derphi0*alpha0);
        b = b / factor;
        // compute the minimizer of the cubic interpolant
        double alpha2 = (-b + sqrt(b*b - 3*a*derphi0)) / (3*a);
        double phi_a2 = phi(n, xk, pk, alpha2);
        // check the break condition
        if (phi_a2 <= phi0 + c1*alpha2*derphi0){
            return std::make_pair(alpha2, phi_a2);
        }
        // update alpha0, alpha1 and alpha2
        if ((alpha1 - alpha2) > alpha1 / 2.0 || (1 - alpha2 / alpha1) < 0.96){
            return std::make_pair(alpha1, phi_a1);
        }
        alpha0 = alpha1;
        alpha1 = alpha2;
        phi_a0 = phi_a1;
        phi_a1 = phi_a2;
    }
    // failed to find a suitable step size
    return std::make_pair(0, phi_a1);
}
// 6. Method line_search_armijo
std::pair<double, double> line_search_armijo( int n, std::vector<double> xk, std::vector<double> pk, std::vector<double> gfk, double old_fval){  
    // get the sum of the product of pk and gfk
    double derphi0 = 0;
    for(int i=0; i<n*n; i++){
        derphi0 += pk[i] * gfk[i];
    }
    std::pair<double, double> pair = scalar_search_armijo(n, xk, pk, old_fval, derphi0);
    return pair;
}

// 7. Method ot_square_regularization
std::vector<double> ot_square_regularization( std::vector<double> a,  std::vector<double> b,  std::vector<double> M, double reg){
    // define G as a.T * b
    int n = a.size();
    std::vector<double> G((int)(n*n));
    for(int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            G[i*n+j] = a[i]*b[j];
        }
    }
    double cost_G = cost(M, G, n, reg);
    double old_cost_G;
    double bound;
    bool stop = false;
    int it = 0;

    while (!stop){
        it += 1;
        old_cost_G = cost_G;
        // compute M + reg*G
        bound = 0;
        std::vector<double> Mi((int)(n*n));
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                Mi[i*n+j] = M[i*n+j] + reg*G[i*n+j];
                if (Mi[i*n+j] < bound){
                    bound = Mi[i*n+j];
                }
            }
        }
        // set M positive
        if (bound < 0){
            for (int i=0; i<n; i++){
                for (int j=0; j<n; j++){
                    Mi[i*n+j] -= bound;
                }
            }
        }

        // solve linear program
        std::vector<double> Gc = emd(a, b, Mi);

        // compute the difference between G and Gc
        std::vector<double> deltaG = std::vector<double>(n*n);
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                deltaG[i*n+j] = Gc[i*n+j] - G[i*n+j];
            }
        }
        std::pair<double, double> pair = line_search_armijo(n, G, deltaG, Mi, cost_G);

        double alpha = pair.first;
        cost_G = pair.second;

        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                G[i*n+j] += alpha*deltaG[i*n+j];
            }
        }

        if (old_cost_G - cost_G < 1e-6){
            stop = true;
        }
        double abs_delta_cost_G = std::abs(old_cost_G - cost_G);
        double rel_delta_cost_G = abs_delta_cost_G / std::abs(cost_G) + 1e-9;
        if (rel_delta_cost_G < 1e-9){
            stop = true;
        }

        // Simetric matrix: we want to ensure that G is symmetric, so we can average G with its transpose (eliminate numerical artefacts)
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) { // Fixa't en el "i + 1"
                double val = (G[i * n + j] + G[j * n + i]) / 2.0;
                G[i * n + j] = val;
                G[j * n + i] = val;
            }
        }

        
    }
    return G;
}

// Allocating the application pid to X cpu
void set_cpu_affinity_V2(uint32_t cpu, pid_t pid)
{
	// Set CPU affinity
	LOGINF("Set CPU {} to PID {}"_format(cpu, pid));
	cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(cpu, &mask);
	if (sched_setaffinity(pid, sizeof(mask), &mask) == -1 ){
		LOGINF("no PID -> {} a cpu"_format(pid, cpu));
		throw_with_trace(std::runtime_error("Could not set CPU affinity: " + std::string(strerror(errno)) ));
	 }

	int new_cpu = get_cpu_id(pid);
	LOGINF("{} {}"_format(new_cpu,cpu));

}
