# fluxpy/auto_decorator.py
import sys
import types
import importlib
from .decorators import print_doc_on_typeerror

def decorate_all_functions_in_module(module):
    for name, obj in vars(module).items():
        if isinstance(obj, types.FunctionType):
            setattr(module, name, print_doc_on_typeerror(obj))

def auto_decorate():
    current_module = sys.modules[__name__.split('.')[0]]
    for name, obj in vars(current_module).items():
        if isinstance(obj, types.ModuleType) and obj.__name__.startswith('fluxpy'):
            decorate_all_functions_in_module(obj)
            # Recursively apply the decorator to submodules
            for submodule_name in dir(obj):
                submodule = getattr(obj, submodule_name)
                if isinstance(submodule, types.ModuleType) and submodule.__name__.startswith('fluxpy'):
                    decorate_all_functions_in_module(submodule)

# Trigger the auto-decorator when the package is imported
auto_decorate()
