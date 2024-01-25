import numpy as np
from pycalphad import variables as v

# In this strategy, convergence failures (Phase all `''`)  have 0 allowed and (1 - 0 = 1) unallowed
def _get_amount_disallowed_phases(eq_result, allowed_phases):
    """
    Determine whether all conditions in the equilibrium result are feasible based on a target phase fraction tolerance

    This function works generally for any sets of conditions, but it probably only makes sense to be used at fixed compositions.

    Parameters
    ----------
    eq_result : xarray.Dataset
    allowed_phases : list
    tol : float
        Maximum allowable phase fraction of allowed phases
    """
    amnt_allowed_phases = eq_result.NP.where(eq_result.Phase.isin(allowed_phases)).sum(dim='vertex')
    amnt_unallowed = (1 - amnt_allowed_phases)
    return amnt_unallowed

# This compute amount unallowed directly
def get_amount_disallowed_phases(eq_result, allowed_phases):
    """
    Determine whether all conditions in the equilibrium result are feasible based on a target phase fraction tolerance

    This function works generally for any sets of conditions, but it probably only makes sense to be used at fixed compositions.

    Parameters
    ----------
    eq_result : xarray.Dataset
    allowed_phases : list
    tol : float
        Maximum allowable phase fraction of allowed phases
    """
    amnt_unallowed = eq_result.NP.where(~eq_result.Phase.isin(allowed_phases)).sum(dim='vertex')
    return amnt_unallowed


# Strategy for this one is to do (1 - allowed phases)
def _get_amount_disallowed_phases_scheil(solidification_result, allowed_phases):

    """
    Determine whether all conditions in the equilibrium result are feasible based on a target phase fraction tolerance

    This function works generally for any sets of conditions, but it probably only makes sense to be used at fixed compositions.

    Parameters
    ----------
    eq_result : xarray.Dataset
    allowed_phases : list
    tol : float
        Maximum allowable phase fraction of allowed phases
    """
    amnt_allowed_phases = sum(solidification_result.cum_phase_amounts.get(phase_name, [0.0])[-1] for phase_name in allowed_phases)
    amnt_unallowed = (1 - amnt_allowed_phases)
    return amnt_unallowed

# The strategy for this one is to count up disallowed phase amounts explictly. That assumes that Scheil would succeeds if it fails completely (usually happens for narrow liquidus)
def get_amount_disallowed_phases_scheil(solidification_result, allowed_phases):
    """
    Determine whether all conditions in the equilibrium result are feasible based on a target phase fraction tolerance

    This function works generally for any sets of conditions, but it probably only makes sense to be used at fixed compositions.

    Parameters
    ----------
    eq_result : xarray.Dataset
    allowed_phases : list
    tol : float
        Maximum allowable phase fraction of allowed phases
    """
    amnt_unallowed = sum(val[-1] for phase_name, val in solidification_result.cum_phase_amounts.items() if phase_name not in allowed_phases)
    return amnt_unallowed


def _build_mass_balanced_grid(num_indep_comps, ngridpts, adjustment=1e-4):
    """Build list of compositions on [0, 1] for each independent component, obeying mass balance"""
    # We subtract twice the adjustment to ensure that we can get near 1.0 without hitting exactly 1.0 (edge of composition space issues)
    grid = np.linspace(0+adjustment, 1-2*adjustment, ngridpts)
    bcast_grid_comps = np.meshgrid(*[grid for _ in range(num_indep_comps)])
    # Compositions where the mass of the independent components is less than unity
    valid_pts = np.nonzero(np.sum(bcast_grid_comps, axis=0) < 1)
    # Index into each independent component's grid to extract out the valid compositions
    # This step converts each grid into 1D instead of len(indep_comps) dimensions.
    grid_comps = [bcast_grid[valid_pts] for bcast_grid in bcast_grid_comps]
    return grid_comps

def _build_composition_list(indep_comps, grid_compositions):
    """
    indep_comps : List[str]
    grid_compositions : ArrayLike
        Should have shape (len(indep_comps), N) with N grid points
    """
    grid_list = np.array(grid_compositions).T.tolist()
    indep_comp_vars = [v.X(ic) for ic in indep_comps]
    return [dict(zip(indep_comp_vars, grid_comps)) for grid_comps in grid_list]
