/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowGeneralisedConcentrations.h"

#include "Conversion.h"

template<>
InputParameters validParams<PorousFlowGeneralisedConcentrations>()
{
  InputParameters params = validParams<PorousFlowMaterialVectorBase>();
  params.addRequiredCoupledVar("generalised_conc_vars", "List of variables that represent the generalised concentrations.  Format is 'f_ph0^c0 f_ph0^c1 f_ph0^c2 ... f_ph0^cN f_ph1^c0 f_ph1^c1 fph1^c2 ... fph1^cN ... fphP^c0 f_phP^c1 fphP^c2 ... fphP^cN' where N=num_primary and P=num_phases. ");
  params.addClassDescription("This Material forms a std::vector<std::vector ...> of generalised concentrations out of the individual generalised concentrations");
  return params;
}

PorousFlowGeneralisedConcentrations::PorousFlowGeneralisedConcentrations(const InputParameters & parameters) :
    PorousFlowMaterialVectorBase(parameters),

    _gen_conc(declareProperty<std::vector<std::vector<Real> > >("PorousFlow_gen_conc")),
    _gen_conc_old(declarePropertyOld<std::vector<std::vector<Real> > >("PorousFlow_gen_conc")),
    _grad_gen_conc(declareProperty<std::vector<std::vector<RealGradient> > >("PorousFlow_grad_gen_conc")),
    _dgen_conc_dvar(declareProperty<std::vector<std::vector<std::vector<Real> > > >("dPorousFlow_gen_conc_dvar")),

    _num_passed_gc_vars(coupledComponents("generalised_conc_vars"))
{  if (_num_passed_gc_vars != _num_phases*_num_primary)
    mooseError("PorousFlowGeneralisedConcentrations: The number of gen_conc_vars is " << _num_passed_gc_vars << " which must be equal to the Dictator's num_phases (" << _num_phases << ") multiplied by num_primary (" << _num_primary << ")");

  _gc_vars_num.resize(_num_passed_gc_vars);
  _gc_vars.resize(_num_passed_gc_vars);
  _grad_gc_vars.resize(_num_passed_gc_vars);
  for (unsigned i = 0; i < _num_passed_gc_vars; ++i)
  {
    _gc_vars_num[i] = coupled("generalised_conc_vars", i);
    _gc_vars[i] = &coupledNodalValue("generalised_conc_vars", i);
    _grad_gc_vars[i] = &coupledGradient("generalised_conc_vars", i);
  }
}

void
PorousFlowGeneralisedConcentrations::initQpStatefulProperties()
{
  _gen_conc[_qp].resize(_num_phases);
  _gen_conc_old[_qp].resize(_num_phases);
  _grad_gen_conc[_qp].resize(_num_phases);
  _dgen_conc_dvar[_qp].resize(_num_phases);
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _gen_conc[_qp][ph].resize(_num_primary);
    _gen_conc_old[_qp][ph].resize(_num_primary);
    _grad_gen_conc[_qp][ph].resize(_num_primary);
    _dgen_conc_dvar[_qp][ph].resize(_num_primary);
    for (unsigned int comp = 0; comp < _num_primary; ++comp)
      _dgen_conc_dvar[_qp][ph][comp].assign(_num_var, 0.0);
  }

  // the derivative matrix is fixed for all time
  // so it can be built here instead of in computeQpProperties
  /**
   * Note that in certain kernels (such as the dispersion kernel)
   * we assume that dgen_conc_dvar = d(grad(gen_conv))/(d(grad(var)))
   */
  unsigned int i = 0;
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    for (unsigned int comp = 0; comp < _num_primary; ++comp)
    {
      if (_dictator.isPorousFlowVariable(_gc_vars_num[i]))
      {
        // _gc_vars[i] is a PorousFlow variable
        const unsigned int pf_var_num = _dictator.porousFlowVariableNum(_gc_vars_num[i]);
        _dgen_conc_dvar[_qp][ph][comp][pf_var_num] = 1.0;
      }
      i++;
    }
  }

  build_gen_conc(_qp);
}

void
PorousFlowGeneralisedConcentrations::computeQpProperties()
{
  build_gen_conc(_qp);
}

void
PorousFlowGeneralisedConcentrations::build_gen_conc(unsigned int qp)
{
  unsigned int i = 0;
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {

    for (unsigned int comp = 0; comp < _num_primary ; ++comp)
    {
      _gen_conc[qp][ph][comp] = (*_gc_vars[i])[_node_number[qp]];
      _grad_gen_conc[qp][ph][comp] = (*_grad_gc_vars[i])[qp];
      i++;
    }
  }
}
