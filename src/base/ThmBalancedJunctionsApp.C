#include "ThmBalancedJunctionsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ThmBalancedJunctionsApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

ThmBalancedJunctionsApp::ThmBalancedJunctionsApp(const InputParameters & parameters)
  : MooseApp(parameters)
{
  ThmBalancedJunctionsApp::registerAll(_factory, _action_factory, _syntax);
}

ThmBalancedJunctionsApp::~ThmBalancedJunctionsApp() {}

void
ThmBalancedJunctionsApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<ThmBalancedJunctionsApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"ThmBalancedJunctionsApp"});
  Registry::registerActionsTo(af, {"ThmBalancedJunctionsApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ThmBalancedJunctionsApp::registerApps()
{
  registerApp(ThmBalancedJunctionsApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ThmBalancedJunctionsApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ThmBalancedJunctionsApp::registerAll(f, af, s);
}
extern "C" void
ThmBalancedJunctionsApp__registerApps()
{
  ThmBalancedJunctionsApp::registerApps();
}
