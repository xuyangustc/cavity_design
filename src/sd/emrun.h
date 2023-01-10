/*
 * emrun.h
 *
 *  Created on: 2018年7月5日
 *      Author: hyliu
 */

#ifndef EMRUN_H_
#define EMRUN_H_
#include "sd/forcefield.h"
#include "gsl/gsl_multimin.h"
#include "Eigen/Core"
//#include "LBFGS.h"
#include "LBFGS_addstop.h"
namespace NSPsd {
class EMRun;
void enededx(const gsl_vector *x,void *emrun,double *ene,gsl_vector *dedx);
double ene(const gsl_vector *x, void *emrun);
void dedx(const gsl_vector *x,void *emrun,gsl_vector *dedx);

class EneFunctor {
public:
	EneFunctor(int n, EMRun *emrun, int si, double ed, int ms, const std::vector<bool> &fce, const std::vector<double> &ro):
		n_(n),emrun_(emrun),step_interval(si),energy_diff(ed),max_step(ms),forceoff_(fce),rmsd_ori(ro){;}
	double operator()(const Eigen::VectorXd &x, Eigen::VectorXd& grad);

	void add_rmsd();
	bool min() {return min_;}
	std::vector<double> energy() {return min_energies;}
	std::vector<double> crd() {return min_crd;}
	Eigen::VectorXd minx() {return min_x;}
	Eigen::VectorXd mingrad() {return min_grad;}
	//void resetemrunene();
	int nowstep() {return now_step;}
	std::vector<double> enecurve() {return ene_curve;}
	std::vector<double> rmsdmc() {return rmsd_mc;}
	std::vector<double> rmsdall() {return rmsd_all;}
private:
	EMRun *emrun_;
	int n_;
	int max_step{-1};
	const std::vector<bool> forceoff_;
	const std::vector<double> rmsd_ori;

	bool min_{false};
	int now_step{-1};
	int step_interval{0};
	double energy_diff{0.0};

	std::vector<double> records_; //minimum energy

	int min_step{0};
	std::vector<double> min_energies;
	std::vector<double> min_crd;
	Eigen::VectorXd min_x;
	Eigen::VectorXd min_grad;

	std::vector<double> ene_curve;
	std::vector<double> rmsd_mc;
	std::vector<double> rmsd_all;
	//bool rmsdrecord{true};
};
class GenChain;
class EMRun {
public:
	struct Param{
		double gsl_initstepsize{0.1};
		double gsl_linesearchtol{0.1};
		double lbfsg_epsilon{1e-6};
	};
	Param param;
	enum {BFGS2,BFGS,CONJUGATE_PR,CONJUGATE_FR,SDESCENT,LBFGS};
	EMRun(){;}
	EMRun(const GenChain *genchain,const ForceField *ff,
			std::vector<bool> fixedatom=std::vector<bool>(), std::vector<double> natcrd=std::vector<double>());
	bool run(int algorithm, int maxsteps, std::vector<double> *outcrd, double *ene, int step_interval,
			double energy_diff, std::vector<double> *gradient=nullptr){
		if(algorithm==LBFGS)
			return run_lbfgs(maxsteps,outcrd,ene,step_interval,energy_diff,gradient);
		else
			return run_gsl(algorithm,maxsteps,outcrd,ene,gradient);
	}
	void writepdb(std::ostream &ofs);
	bool run_gsl(int algorithm,int maxsteps,std::vector<double> *outcrd,double *ene,
			std::vector<double> *gradient);
	bool run_lbfgs(int maxsteps,std::vector<double> *outcrd,double *ene,int step_interval,
			double energy_diff,std::vector<double> *gradient=nullptr);
	const ForceField *forcefield() const {return ff_;}
	std::vector<double> &crdbck() {return crdbck_;}
	std::vector<bool> &atomfixed(){return atomfixed_;}
	std::vector<double> &potenergies() {return potenergies_;}
	std::vector<double> &forces() {return forces_;}
	ActiveSelections &acts(){ return *acts_;}

	std::vector<double> enecurve() {return ene_curve;}
	//std::vector<bool> & atomfixed() {return atomfixed_;}
	std::vector<double> &nativecrd() {return native_crd;}
	std::vector<double> rmsdmc() {return rmsd_mc;}
	std::vector<double> rmsdall() {return rmsd_all;}
protected:
	const ForceField* ff_;
	const GenChain *genchain_;
	std::vector<bool> atomfixed_;
	std::vector<double> crdbck_;
	std::vector<double> potenergies_;
	std::vector<double> forces_;
	ActiveSelections *acts_;
	int ndof_;
	void set_gsl_vector(const std::vector<double> &source, gsl_vector *dest);
	void copy_gsl_vector(const gsl_vector *source, std::vector<double> &dest);

	std::vector<double> ene_curve;
	std::vector<double> native_crd;
	std::vector<double> rmsd_mc;
	std::vector<double> rmsd_all;
};

}



#endif /* EMRUN_H_ */
