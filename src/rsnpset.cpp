#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

MatrixXd getU(const VectorXd& Y, const VectorXd& delta, const MatrixXd& G,const std::string& score) {
    if(score=="cox") {
        const int n = Y.rows();
        MatrixXd YY   = Y.rowwise().replicate(n);
        MatrixXd matY = (YY.transpose().array()>=YY.array()).select(MatrixXd::Ones(n,n),0);
        MatrixXd matA = matY*G;

        VectorXd vecYsum = matY.rowwise().sum();
        MatrixXd U = G-(vecYsum.cwiseInverse().asDiagonal()*matA);
        return delta.asDiagonal()*U;
    }
    else { // else if(score=="binomial") { // else if(score=="gaussian") {
        VectorXd Ytilde = Y.array() - Y.mean();
        MatrixXd Gtilde = G-G.colwise().mean().colwise().replicate(G.rows());
        return Ytilde.asDiagonal()*Gtilde;
    }
}

MatrixXd submatcols(const MatrixXd& m, const ArrayXi& idx ) {
    MatrixXd subm(m.rows(), idx.size());
    for(int i=0;i<idx.size();i++) {
        subm.col(i)=m.col(idx[i]);
    }
    return subm;
}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

RcppExport SEXP rsnpsetRcpp (const SEXP Y_, const SEXP delta_, const SEXP G_, const SEXP geneIdxSets_, 
                             const SEXP permutation_, const SEXP score_, const SEXP vMethod_, 
                             const SEXP pinvCheck_, const SEXP tolerance_) {
    std::string score = as<std::string>(score_); 
    std::string vMethod   = as<std::string>(vMethod_); 
    double tolerance = as<double>(tolerance_); 
    bool permutation = as<bool>(permutation_); 
    bool pinvCheck   = as<bool>(pinvCheck_); 
    const Map<VectorXd> Y    (as<Map<VectorXd> >(Y_));
    const Map<VectorXd> delta(as<Map<VectorXd> >(delta_));
    const Map<MatrixXd> G    (as<Map<MatrixXd> >(G_));

    MatrixXd U;
    if(permutation==true) {
        RNGScope scope;
        int n = Y.size();
        PermutationMatrix<Dynamic,Dynamic> perm(n);
        perm.setIdentity();
        std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(),randWrapper);
        VectorXd Yperm     = perm * Y;
        VectorXd deltaPerm = perm * delta;
        U = getU(Yperm,deltaPerm,G,score);
    }
    else{
        U = getU(Y,delta,G,score);
    }
    
    List geneIdxSets(geneIdxSets_);
    const int geneSetNum = geneIdxSets.size();
    
    ArrayXd __W(geneSetNum);
    ArrayXi _rk(geneSetNum);
    ArrayXd _d0(geneSetNum);
    ArrayXd _d1(geneSetNum);
    ArrayXd _d2(geneSetNum);
    ArrayXd _d3(geneSetNum);
    ArrayXd _d4(geneSetNum);
    
    double muhat = Y.sum()/Y.size();
    double muhatpq = muhat*(1-muhat);
    MatrixXd centeredG = G - (MatrixXd::Ones(G.rows(),1) * G.colwise().mean() ) ;
    
    for(int i=0;i<geneSetNum;i++) {     
        ArrayXi geneSet = (as<Map<ArrayXi> >(geneIdxSets[i]))-1; // idx diff betw R & C
        
        MatrixXd Usub = submatcols(U ,geneSet);
        MatrixXd V;
        if(vMethod=="empirical") {
            V = Usub.transpose() * Usub; 
        }
        else if(score=="binomial" && vMethod=="asymptotic") {
            MatrixXd Gsub = submatcols(centeredG, geneSet);
            MatrixXd GtG  = Gsub.transpose()*Gsub;
            V = muhatpq*GtG;
        }
        
        SelfAdjointEigenSolver<MatrixXd> sdecomp(V);
        ArrayXd evals=sdecomp.eigenvalues();
        MatrixXd Q=sdecomp.eigenvectors();

        double maxeval=evals.maxCoeff();
        double tol=maxeval*tolerance;
        const int rk((evals>tol).count());
        const ArrayXd evalsp=evals.tail(rk);
        MatrixXd Qp=Q.rightCols(rk);
        MatrixXd Dp(evalsp.matrix().asDiagonal());
        MatrixXd Dpinv((1/evalsp).matrix().asDiagonal());
        MatrixXd VPM=Qp*Dpinv*Qp.transpose();
        
        MatrixXd U0 = Usub.colwise().sum().transpose();
        MatrixXd Wmatrix= U0.transpose()* VPM *U0;
        __W[i]=(Wmatrix)(0,0);
        _rk[i]=rk;
        
        if(pinvCheck) {
            _d0[i]=(V-Qp*Dp*Qp.transpose()).lpNorm<Infinity>(); 
            _d1[i]=(V*VPM*V-V).lpNorm<Infinity>(); 
            _d2[i]=(VPM*V*VPM-VPM).lpNorm<Infinity>(); 
            _d3[i]=((V*VPM).transpose()-V*VPM).lpNorm<Infinity>();
            _d4[i]=((VPM*V).transpose()-VPM*V).lpNorm<Infinity>();
        }
    }
    if(pinvCheck) {
        return DataFrame::create( Named("W"   )=wrap(__W), 
                                           Named("rank")=wrap(_rk),
                                           Named("d0"  )=wrap(_d0),
                                           Named("d1"  )=wrap(_d1),
                                           Named("d2"  )=wrap(_d2),
                                           Named("d3"  )=wrap(_d3),
                                           Named("d4"  )=wrap(_d4) );
    }
    else {
        return DataFrame::create( Named("W"   )=wrap(__W), 
                                           Named("rank")=wrap(_rk) );
    }
}
