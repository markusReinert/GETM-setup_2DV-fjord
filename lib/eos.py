"""Functions related to the equation of state of sea water.

by Markus Reinert, July 2021
"""

import numpy as np


def rho(SP, PT, p):
    """Density of sea water by Jackett et al. (2006).

    This code is taken from GETM and converted to Python.

    Input parameters:
      SP: (practical/in-situ) salinity (psu)
      PT: potential temperature        (deg C, ITS-90)
      p : gauge pressure               (dbar)
          (absolute pressure - 10.1325 dbar)

    Output value: in-situ density (kg/m³)
    
    Check value: rho(20, 20, 1000) = 1017.728868019642

    Comments from the Fortran code:
    based on DRJ on 10/12/03
    REVISION HISTORY:
    AS 2009 based on code provided by Jackett 2005
    See the log for the module
    """
    PT2 = PT ** 2

    anum = (
        9.9984085444849347e02
        + PT
        * (
            7.3471625860981584e00
            + PT * (-5.3211231792841769e-02 + PT * 3.6492439109814549e-04)
        )
        + SP
        * (
            2.5880571023991390e00
            - PT * 6.7168282786692355e-03
            + SP * 1.9203202055760151e-03
        )
    )

    aden = (
        1
        + PT
        * (
            7.2815210113327091e-03
            + PT
            * (
                -4.4787265461983921e-05
                + PT * (3.3851002965802430e-07 + PT * 1.3651202389758572e-10)
            )
        )
        + SP
        * (
            1.7632126669040377e-03
            - PT * (8.8066583251206474e-06 + PT2 * 1.8832689434804897e-10)
            + np.sqrt(SP) * (5.7463776745432097e-06 + PT2 * 1.4716275472242334e-09)
        )
    )

    density = anum / aden

    if p != 0:
        pPT = p * PT
        anum = anum + p * (
            1.1798263740430364e-02
            + PT2 * 9.8920219266399117e-08
            + SP * 4.6996642771754730e-06
            - p * (2.5862187075154352e-08 + PT2 * 3.2921414007960662e-12)
        )
        aden = aden + p * (
            6.7103246285651894e-06
            - pPT * (PT2 * 2.4461698007024582e-17 + p * 9.1534417604289062e-18)
        )
        density = anum / aden

    return density


lambda_1 = -5.73e-2
lambda_2 = 8.32e-2
lambda_3 = 7.61e-4
def T_freezing(S, z):
    """Linearized freezing temperature of sea water by Jenkins (2011)."""
    return lambda_1 * S + lambda_2 + lambda_3 * z
