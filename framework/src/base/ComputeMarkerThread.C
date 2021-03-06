/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

// MOOSE includes
#include "ComputeMarkerThread.h"
#include "AuxiliarySystem.h"
#include "Problem.h"
#include "FEProblem.h"
#include "Marker.h"

// libmesh includes
#include "libmesh/threads.h"

ComputeMarkerThread::ComputeMarkerThread(FEProblem & fe_problem, AuxiliarySystem & sys) :
    ThreadedElementLoop<ConstElemRange>(fe_problem, sys),
    _fe_problem(fe_problem),
    _aux_sys(sys),
    _marker_whs(_fe_problem.getMarkerWarehouse())
{
}

// Splitting Constructor
ComputeMarkerThread::ComputeMarkerThread(ComputeMarkerThread & x, Threads::split split) :
    ThreadedElementLoop<ConstElemRange>(x, split),
    _fe_problem(x._fe_problem),
    _aux_sys(x._aux_sys),
    _marker_whs(x._marker_whs)
{
}

ComputeMarkerThread::~ComputeMarkerThread()
{
}

void
ComputeMarkerThread::subdomainChanged()
{
  _fe_problem.subdomainSetup(_subdomain, _tid);
  _marker_whs.subdomainSetup(_tid);

  std::set<MooseVariable *> needed_moose_vars;
  _marker_whs.updateVariableDependency(needed_moose_vars, _tid);

  for (std::map<std::string, MooseVariable *>::iterator it = _aux_sys._elem_vars[_tid].begin(); it != _aux_sys._elem_vars[_tid].end(); ++it)
  {
    MooseVariable * var = it->second;
    var->prepareAux();
  }

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
}

void
ComputeMarkerThread::onElement(const Elem *elem)
{
  _fe_problem.prepare(elem, _tid);
  _fe_problem.reinitElem(elem, _tid);
  _fe_problem.reinitMaterials(_subdomain, _tid);

  if (_marker_whs.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<MooseSharedPointer<Marker> > & markers = _marker_whs.getActiveBlockObjects(_subdomain, _tid);
    for (std::vector<MooseSharedPointer<Marker> >::const_iterator it = markers.begin(); it != markers.end(); ++it)
      (*it)->computeMarker();
  }

  _fe_problem.swapBackMaterials(_tid);

  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (std::map<std::string, MooseVariable *>::iterator it = _aux_sys._elem_vars[_tid].begin(); it != _aux_sys._elem_vars[_tid].end(); ++it)
    {
      MooseVariable * var = it->second;
      var->insert(_aux_sys.solution());
    }
  }
}

void
ComputeMarkerThread::onBoundary(const Elem * /*elem*/, unsigned int /*side*/, BoundaryID /*bnd_id*/)
{
}

void
ComputeMarkerThread::onInternalSide(const Elem * /*elem*/, unsigned int /*side*/)
{
}

void
ComputeMarkerThread::postElement(const Elem * /*elem*/)
{
}

void
ComputeMarkerThread::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
}

void
ComputeMarkerThread::join(const ComputeMarkerThread & /*y*/)
{
}
