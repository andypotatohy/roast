"""Makes ``tensorflow.keras`` resolve to legacy Keras 2 before TF is imported.

The bundled ``lib/multiaxial/*.h5`` models were saved under Keras 2 (TF
2.11, per ``lib/multiaxial/multiaxialEnvLinux.yml``). Loading them under
Keras 3 (the default since TF 2.16) fails with e.g.::

    TypeError: Unrecognized keyword arguments passed to Conv2DTranspose: {'groups': 1}

because Keras 3's layer configs dropped some Keras-2-only constructor
kwargs that got serialized into these older HDF5 files. The documented fix
is the ``tf-keras`` compatibility package plus the ``TF_USE_LEGACY_KERAS``
environment variable, which must be set *before* ``tensorflow`` is first
imported anywhere in the process (Keras resolves this at import time, not
per-call). Verified against the actual bundled ``sagittal_model.h5``.

Import this module before importing ``tensorflow`` anywhere else in the
process (``roast_py.segmentation.multiaxial`` does this for you).
"""

from __future__ import annotations

import os
import sys


def ensure_legacy_keras() -> None:
    if "tensorflow" in sys.modules:
        # Keras has already resolved; setting the env var now would be a no-op.
        # Only warn if it resolved to the incompatible Keras 3.
        import tensorflow as tf

        if not type(tf.keras.models.Model).__module__.startswith("tf_keras"):
            raise RuntimeError(
                "tensorflow was already imported with Keras 3 active, which "
                "cannot load the legacy Keras-2 .h5 models bundled under "
                "lib/multiaxial/. Set TF_USE_LEGACY_KERAS=1 in the "
                "environment before starting Python (or before any other "
                "code imports tensorflow), and `pip install tf-keras`."
            )
        return

    os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")
    try:
        import tf_keras  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "roast_py's multiaxial segmentation needs the legacy Keras 2 "
            "runtime to load the bundled .h5 models. Install it with "
            "`pip install tf-keras` (or `pip install roast_py[multiaxial]`)."
        ) from e
