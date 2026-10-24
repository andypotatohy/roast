"""ROAST's SPM routines through SPM's official Python interface, spm-python.

    python -m roast_py.spm status                    # what's installed
    python -m roast_py.spm install --accept-license  # Python packages + MATLAB Runtime
    python -m roast_py.spm check                     # start SPM once and print its version

roast() calls these itself when it needs SPM (multiaxial off, a T2 image,
or resampling); see runtime.py for how the MATLAB Runtime is found or
installed, and routines.py for the SPM steps.
"""
