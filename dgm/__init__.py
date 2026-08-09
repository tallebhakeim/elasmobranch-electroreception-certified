"""dgm (frozen subset): the discrete geometric method with certified complementary bounds.

This is a snapshot of the modules needed to reproduce the results of

    An error-bounded field model of elasmobranch electroreception: canal funnelling,
    frequency tuning and the effect of body plan in a hammerhead shark and a manta ray

It is not the full engine, and it is not maintained here; it is frozen so that the article
stays reproducible independently of later development.

Method: the potential lives on the primal mesh and the current flux on its circumcentric
dual. The primal energy over-estimates the dissipation and the RT0 complementary flux
under-estimates it, so the exact value is enclosed. See the article's supplementary
information, section S3.
"""
__version__ = "1.0.0-article"

from .mesh2d import TriMesh
from .primal import EPS0, assemble_primal, solve_dirichlet, energy
from .mesh3d import TetMesh
from .primal3d import assemble_primal_3d, element_gradients
from .mixed2d import energy_lower_2d

__all__ = ["TriMesh", "TetMesh", "EPS0", "assemble_primal", "solve_dirichlet",
           "energy", "assemble_primal_3d", "element_gradients", "energy_lower_2d"]
