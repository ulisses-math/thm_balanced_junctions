//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ThmBalancedJunctionsTestApp.h"
#include "ThmBalancedJunctionsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
ThmBalancedJunctionsTestApp::validParams()
{
  InputParameters params = ThmBalancedJunctionsApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

ThmBalancedJunctionsTestApp::ThmBalancedJunctionsTestApp(const InputParameters & parameters)
  : MooseApp(parameters)
{
  ThmBalancedJunctionsTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ThmBalancedJunctionsTestApp::~ThmBalancedJunctionsTestApp() {}

void
ThmBalancedJunctionsTestApp::registerAll(Factory & f,
                                         ActionFactory & af,
                                         Syntax & s,
                                         bool use_test_objs)
{
  ThmBalancedJunctionsApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ThmBalancedJunctionsTestApp"});
    Registry::registerActionsTo(af, {"ThmBalancedJunctionsTestApp"});
  }
}

void
ThmBalancedJunctionsTestApp::registerApps()
{
  registerApp(ThmBalancedJunctionsApp);
  registerApp(ThmBalancedJunctionsTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ThmBalancedJunctionsTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ThmBalancedJunctionsTestApp::registerAll(f, af, s);
}
extern "C" void
ThmBalancedJunctionsTestApp__registerApps()
{
  ThmBalancedJunctionsTestApp::registerApps();
}
