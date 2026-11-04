"""Printed final equation, independently expressed without the symbolic engine.

Inputs are a DataFrame with the unsuffixed names returned by chiral_sr.data.load.
Concentration: mg/mL; Mn: g/mol; S, B, R, U, Q, E: mol/kg;
H: solvent H-bond acceptor count; F: log10 backbone fragment solubility in mol/L.
Coefficients use the numeric values in exactly these units. No standardization.
"""
import numpy as np
from scipy.special import expit

B0 = -0.611886980460750
B1 = 0.000920944987798687
B2 = 0.713484255385024
B3 = 0.0000799826488026122
B4 = -0.00206753621927866
B5 = 0.0000000502316279772739


def logit(X):
    c = np.asarray(X['Concentration_mg_mL'], float)
    M = np.asarray(X['MW_Mn'], float)
    S = np.asarray(X['SC_C_atoms'], float)
    B = np.asarray(X['BB_C_atoms'], float)
    R = np.asarray(X['SC_n_branch_points'], float)
    U = np.asarray(X['SC_n_units'], float)
    H = np.asarray(X['Solvent_NumHAcceptors'], float)
    F = np.asarray(X['Backbone_fastsolv'], float)
    Q = np.asarray(X['BB_conjugated_units'], float)
    E = np.asarray(X['BB_hetero_atoms'], float)
    if not all(np.isfinite(v).all() for v in [c, M, S, B, R, U, H, F, Q, E]):
        raise ValueError('All descriptors must be finite; no implicit imputation.')
    if np.any(Q <= 0) or np.any(E <= 0) or np.any(M <= 0) or np.any(c < 0):
        raise ValueError('Require positive Q, E, M and nonnegative c.')
    return B0 + B1*(S-B)*c*R + B2*H*U*F + B3*c*c*F/Q + B4*(S-B)*U*B + B5*c*M*B/E


def predict(X):
    return expit(logit(X))
