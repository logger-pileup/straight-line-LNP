import sys
sys.path.insert(1, './lattice-estimator')
from estimator import *
from estimator.schemes import LWEParameters, SISParameters
from estimator.nd import Uniform

def sis_estimate(n, d, logq, B):
    """
    Wrapper for SIS hardness estimator.

    Parameters
    ----------
    n : length of SIS solution.
    d : dimension of Rq.
    logq : number of bits of modulus q.
    B : L2 norm bound on SIS solution.

    Returns
    -------
    The log (number of bits) of the smallest estimates number of operations
    of any of the attacks on SIS with the given parameters.

    If no attacks were provided, then it returns oo.
    """
    est = SIS.estimate(SISParameters(n*d, 2^logq, B))
    return log(min([L["rop"] for L in est.values()], default=oo), 2).n()


def lwe_estimate(n, d, logq, nu, m=oo):
    """
    Wrapper for LWE hardness estimator.

    Parameters
    ----------
    n : dimensiom of LWE sample.
    d : dimension of Rq.
    logq : number of bits of modulus q.
    nu : maximum size of coefficients of LWE secret and error.
    m : number of samples.

    Returns
    -------
    The log (number of bits) of the smallest estimates number of operations
    of any of the attacks on LWE with the given parameters.

    If no attacks were provided, then it returns oo.
    """
    est = LWE.estimate(LWEParameters(n*d, 2^logq, Uniform(-nu,nu), 
                                     Uniform(-nu,nu), m=m*d), jobs=8)
    return log(min([L["rop"] for L in est.values()], default=oo), 2).n()

def print_dashes(text="", total=100):
        """
        Print a dashed line with some text in it.
        """
        start = (total-len(text))//2
        end = total - len(text) - start
        print(start*"-" + text + end*"-")
