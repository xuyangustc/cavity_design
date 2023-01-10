/*
 * emrun.cpp
 *
 *  Created on: 2018年7月5日
 *      Author: hyliu
 */
#include "sd/emrun.h"
#include "sd/genchain.h"
#include <cmath>
using namespace NSPsd;
void NSPsd::enededx(const gsl_vector *x, void *emrun, double *ene,
		gsl_vector *dedx) {
	const ForceField &ff = *(((EMRun *) emrun)->forcefield());
	std::vector<double> &crd = ((EMRun *) emrun)->crdbck();
	std::vector<bool> & fixed = ((EMRun *) emrun)->atomfixed();
	std::vector<double> &potenergies = ((EMRun*) emrun)->potenergies();
	int natoms = ff.natoms();
	int idx = 0;
	for (int i = 0; i < natoms; ++i) {
		if (fixed[i])
			continue;
		for (int m = 0; m < 3; ++m)
			crd[3 * i + m] = gsl_vector_get(x, idx++);
	}
	NeighborList nbl(crd, ff, &fixed);
	double maxsteric = StericAtom::DEDRMAX;
	StericAtom::DEDRMAX = 1.e32;
	std::vector<double> & forces =((EMRun *) emrun)->forces();
	forces= ff.forces(crd, nbl, &potenergies, ((EMRun *) emrun)->acts());
	StericAtom::DEDRMAX = maxsteric;
	idx = 0;
	for (int i = 0; i < natoms; ++i) {
		if (fixed[i])
			continue;
		for (int m = 0; m < 3; ++m)
			gsl_vector_set(dedx, idx++, -forces[3 * i + m]);
	}
	*ene = potenergies[ForceField::ETOT];
}
void NSPsd::dedx(const gsl_vector *x, void *emrun, gsl_vector *dedx) {
	double ene;
	enededx(x, emrun, &ene, dedx);
}
double NSPsd::ene(const gsl_vector *x, void *emrun) {
	double ene;
	gsl_vector *dedx = gsl_vector_alloc(x->size);
	enededx(x, emrun, &ene, dedx);
	gsl_vector_free(dedx);
	return ene;
}
void EneFunctor::add_rmsd() {
	rmsd_mc.push_back(emrun_->forcefield()->getmcrmsd(emrun_->crdbck(),rmsd_ori,forceoff_,10.0));
	rmsd_all.push_back(emrun_->forcefield()->getallrmsd(emrun_->crdbck(),rmsd_ori,forceoff_,10.0));
}
double EneFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
	//static int ncalls{0};
	auto gsl_x = gsl_vector_alloc(n_);
	auto gsl_grad = gsl_vector_alloc(n_);
	for (int i = 0; i < n_; ++i)
		gsl_vector_set(gsl_x, i, x[i]);
	double e;
	enededx(gsl_x, (void *) emrun_, &e, gsl_grad);
	for (int i = 0; i < n_; ++i)
		grad[i] = gsl_vector_get(gsl_grad, i);

	if(++now_step % 20==0){
		double gnorm=0.0;
		for(int i=0;i<n_;++i){
			gnorm +=grad[i]*grad[i];
		}
		gnorm=sqrt(gnorm/(double) n_);
		//std::cout << "Number_ene_calculations=" << now_step
		//		<<"  Current_ene=" << emrun_->potenergies()[ForceField::ETOT]
		//		<<"  RMS_gradient=" <<gnorm	<<std::endl;
		//std::ofstream ofs("outem.pdb");
		//emrun_->writepdb(ofs);
		//ofs.close();
	}
	//++now_step;

	gsl_vector_free(gsl_x);
	gsl_vector_free(gsl_grad);

	if(std::isnan(e)) {
		min_ = true;
		std::cout <<"Optimization:  " <<now_step <<"  " <<min_step <<"  " <<min_energies[0] <<std::endl;
		return e;
	}
	ene_curve.push_back(e);
	if(!rmsd_ori.empty()) add_rmsd();
	if(records_.empty() || e<records_.back()) {
		records_.push_back(e);
		min_step = now_step;
		min_energies = emrun_->potenergies();
		min_crd = emrun_->crdbck();
		min_x = x;
		min_grad = grad;
	} else records_.push_back(records_.back());
	int prev = records_.size() - step_interval;
	if(prev>=0) {
		if(records_[prev]-records_.back()<energy_diff) {
			min_ = true;
			std::cout <<"Optimization:  " <<now_step <<"  " <<min_step <<"  " <<min_energies[0] <<std::endl;
			return e;
		}
	}
	if(now_step>=max_step) {
		min_ = true;
		std::cout <<"Optimization:  " <<now_step <<"  " <<min_step <<"  " <<min_energies[0] <<std::endl;
		return e;
	}
	return e;
}

EMRun::EMRun(const GenChain * genchain,const ForceField* ff, std::vector<bool> fixed, std::vector<double> natcrd) :
		genchain_(genchain),ff_(ff) ,atomfixed_(fixed), native_crd(natcrd){
	acts_=genchain->getactiveselections();
	std::vector<bool> afix=acts_->atomfixed();
	if (atomfixed_.empty()) {
		atomfixed_.assign(ff_->natoms(), false);
	}
	for(int i=0;i<ff_->natoms();++i){
		atomfixed_[i]= atomfixed_[i] || afix[i];
	}
	ndof_ = 3 * ff_->natoms();
	for (auto f : atomfixed_)
			if (f)	ndof_-=3;
}

void EMRun::writepdb(std::ostream &os){
	genchain_->writepdb(crdbck_,os,1.0/A2NM);
}
void EMRun::set_gsl_vector(const std::vector<double> &crd, gsl_vector *x) {
	int idx = 0;
	for (int i = 0; i < crd.size() / 3; ++i) {
		if (atomfixed_[i])
			continue;
		for (int m = 0; m < 3; ++m)
			gsl_vector_set(x, idx++, crd[3 * i + m]);
	}
}
void EMRun::copy_gsl_vector(const gsl_vector *x, std::vector<double> &crd) {
	int idx = 0;
	for (int i = 0; i < crd.size() / 3; ++i) {
		if (atomfixed_[i])
			continue;
		for (int m = 0; m < 3; ++m)
			crd[3 * i + m] = gsl_vector_get(x, idx++);
	}
}

bool EMRun::run_gsl(int algorithm, int maxsteps, std::vector<double> *outcrd, double *outene,
		std::vector<double> *gradient) {
	static std::map<int, const gsl_multimin_fdfminimizer_type *> algorithms { {
			BFGS2, gsl_multimin_fdfminimizer_vector_bfgs2 }, { BFGS,
			gsl_multimin_fdfminimizer_vector_bfgs }, { CONJUGATE_FR,
			gsl_multimin_fdfminimizer_conjugate_fr }, { CONJUGATE_PR,
			gsl_multimin_fdfminimizer_conjugate_pr }, { SDESCENT,
			gsl_multimin_fdfminimizer_steepest_descent } };
	crdbck_ = *outcrd;
	gsl_vector *x = gsl_vector_alloc(ndof_);
	auto minimizer = gsl_multimin_fdfminimizer_alloc(algorithms.at(algorithm),
			ndof_);
	set_gsl_vector(crdbck_, x);
	gsl_multimin_function_fdf fdf;
	fdf.n = ndof_;
	fdf.f = ene;
	fdf.df = dedx;
	fdf.fdf = enededx;
	fdf.params = (void *) this;
	gsl_multimin_fdfminimizer_set(minimizer, &fdf, x, param.gsl_initstepsize,
			param.gsl_linesearchtol);
	int iter = 0;
	int status;
	while (iter < maxsteps) {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(minimizer);
		if (status == GSL_ENOPROG)
			break;
		status = gsl_multimin_test_gradient(minimizer->gradient, 1e-3);
		if (status == GSL_SUCCESS)
			break;
		std::cout << "Iteration " << iter << ": "
				<< gsl_multimin_fdfminimizer_minimum(minimizer) << std::endl;
	}
	copy_gsl_vector(gsl_multimin_fdfminimizer_x(minimizer), crdbck_);
	*outcrd = crdbck_;
	*outene = gsl_multimin_fdfminimizer_minimum(minimizer);
	if (gradient) {
		gradient->assign(3 * ff_->natoms(), 0.0);
		copy_gsl_vector(minimizer->gradient, *gradient);
	}
	gsl_multimin_fdfminimizer_free(minimizer);
	gsl_vector_free(x);
	if (status == GSL_SUCCESS) {
		return true;
	} else {
		return false;
	}
}
bool EMRun::run_lbfgs(int maxsteps, std::vector<double> *outcrd, double *outene,
		int step_interval, double energy_diff, std::vector<double> *gradient) {
	using namespace LBFGSpp;
	crdbck_ = *outcrd;
	LBFGSParam<double> lbfsg_param;
	lbfsg_param.epsilon = 1e-6;
	lbfsg_param.max_iterations = maxsteps;
	LBFGSSolver<double> solver(lbfsg_param);
	EneFunctor fun(ndof_, this, step_interval, energy_diff, maxsteps, atomfixed_, native_crd);
	Eigen::VectorXd x;
	x.resize(ndof_);
	int idx = 0;
	for (int i = 0; i < crdbck_.size() / 3; ++i) {
		if (atomfixed_[i])
			continue;
		for (int m = 0; m < 3; ++m)
			x[idx++] = crdbck_[3 * i + m];
	}
	double fx;
//test
/*	Eigen::VectorXd grad(ndof_),gradtmp(ndof_);
	double e=fun(x,grad);
	for(int i=0;i<ndof_;++i){
		x[i]=x[i]+0.0001;
		double ep=fun(x,gradtmp);
		x[i]=x[i]-0.0002;
		double em=fun(x,gradtmp);
		std::cout <<i<<" "<<(ep-em)/0.0002<<" "<< grad[i] << std::endl;
		x[i] +=0.0001;
	}*/
	//int niter=0;
	//while (niter<maxsteps){
		int niter_p = solver.minimize(fun, x, fx);
    	x = fun.minx();
    	fx = fun.energy()[0];
    	crdbck_ = fun.crd();
    	potenergies_ = fun.energy();
    	ene_curve.clear();
    	ene_curve = fun.enecurve();
    	rmsd_mc=fun.rmsdmc();
    	rmsd_all=fun.rmsdall();
		//niter+=niter_p;
		//std::cout<<"Number of LBFGS iterations: "<< niter <<std::endl;
		//if(niter_p < 200) break;
	//}
	idx = 0;
	for (int i = 0; i < crdbck_.size() / 3; ++i) {
		if (atomfixed_[i])
			continue;
		for (int m = 0; m < 3; ++m)
			crdbck_[3 * i + m] = x[idx++];
	}
	*outcrd = crdbck_;
	*outene = fx;
	if (gradient) {
		Eigen::VectorXd grad;
		grad.resize(ndof_);
		fun(x, grad);
		idx = 0;
		gradient->clear();
		for (int i = 0; i < crdbck_.size() / 3; ++i) {
			if (atomfixed_[i])
				continue;
			for (int m = 0; m < 3; ++m)
				gradient->push_back(grad(idx++));
		}
	}
	return true;
}
