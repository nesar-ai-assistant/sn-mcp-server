"""Science tool package. Each family is a sub-package with its own __all__;
[tool.mcp-server].tool_modules in pyproject.toml lists which modules the
server registers. Kernel functions declared in a family's DISPATCH_KERNELS
must import only numpy/scipy/astropy + their declared pip_deps."""
