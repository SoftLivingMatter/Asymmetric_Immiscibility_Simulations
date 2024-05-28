import sys
import pytest

# mask hard-to-install dependencies
# still need to patch individually

module = type(sys)('gsd')
module.submodule = type(sys)("hoomd")
sys.modules['gsd'] = module
sys.modules["gsd.hoomd"] = module.submodule

module = type(sys)('hoomd')
module.md = type(sys)("md")
module.azplugins = type(sys)("azplugins")
sys.modules['hoomd'] = module
sys.modules["hoomd.md"] = module.md
sys.modules["hoomd.azplugins"] = module.azplugins
