import pandas as pd

def Solo_M1mHpC4H11N(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + mz - 72.081324
    mz_Cl = 34.968853 + mz - 72.081324
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl])

def Solo_M1mHpHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + mz - 44.997654
    mz_Cl = 34.968853 + mz - 44.997654
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl])

def Solo_M1m2HpNapHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + mz - 66.979600
    mz_Cl = 34.968853 + mz - 66.979600
    mz_m2HpNa = 20.97412 + mz - 66.979600
    
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa])


def Solo_M1m2HpNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + mz - 66.979600
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H])

def Solo_M1m2HpK(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + mz - 36.948058
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H])

def Solo_M2mHpC4H11N(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 72.081324)/2
    mz_Cl = 34.968853 + (mz - 72.081324)/2
    mz_m2HpNa = 20.97412 + (mz - 72.081324)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa])

def Solo_M2mHpHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 44.997654)/2
    mz_Cl = 34.968853 + (mz - 44.997654)/2
    mz_m2HpNa = 20.97412 + (mz - 44.997654)/2
    mz_mHpHCOOH = 44.997654 + (mz - 44.997654)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_mHpHCOOH = peaks.between(mz_mHpHCOOH - prec_mass_error, mz_mHpHCOOH + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_mHpHCOOH])

def Solo_M2mH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz + 1.007825)/2
    mz_Cl = 34.968853 + (mz + 1.007825)/2
    mz_m2HpNa = 20.97412 + (mz + 1.007825)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa])

def Solo_M2pCl(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 34.968853)/2
    mz_Cl = 34.968853 + (mz - 34.968853)/2
    mz_m2HpNa = 20.97412 + (mz - 34.968853)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa])

def Solo_M2m2HpNapHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 66.979600)/2
    mz_Cl = 34.968853 + (mz - 66.979600)/2
    mz_m2HpNa = 20.97412 + (mz - 66.979600)/2
    mz_m2HpNapHCOOH = 66.9796 + (mz - 66.979600)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH])

def Solo_M2m2HpNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 20.97412)/2
    mz_Cl = 34.968853 + (mz - 20.97412)/2
    mz_m2HpNa = 20.97412 + (mz - 20.97412)/2
    mz_m2HpNapHCOOH = 66.9796 + (mz - 20.97412)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH])

def Solo_M2m2HpK(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 36.948058)/2
    mz_Cl = 34.968853 + (mz - 36.948058)/2
    mz_m2HpNa = 20.97412 + (mz - 36.948058)/2
    mz_m2HpNapHCOOH = 66.9796 + (mz - 36.948058)/2
    mz_m2HpK = 36.948058 + (mz - 36.948058)/2
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK])

def Solo_M3mH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz + 1.007825)/3
    mz_Cl = 34.968853 + (mz + 1.007825)/3
    mz_m2HpNa = 20.97412 + (mz + 1.007825)/3
    mz_m2HpNapHCOOH = 66.9796 + (mz + 1.007825)/3
    mz_m2HpK = 36.948058 + (mz + 1.007825)/3

    mz_M2mH = -1.007825 + (mz + 1.007825)*(2/3)
    mz_M2pCl = 34.968853 + (mz + 1.007825)*(2/3)
    mz_M2m2HpNa = 20.97412 + (mz + 1.007825)*(2/3)
    mz_M2m2HpNapHCOOH = 66.9796 + (mz + 1.007825)*(2/3)
    mz_M2m2HpK = 36.948058 + (mz + 1.007825)*(2/3)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])

def Solo_M3pCl(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 34.968853)/3
    mz_Cl = 34.968853 + (mz - 34.968853)/3
    mz_m2HpNa = 20.97412 + (mz - 34.968853)/3
    mz_m2HpNapHCOOH = 66.9796 + (mz - 34.968853)/3
    mz_m2HpK = 36.948058 + (mz - 34.968853)/3

    mz_M2mH = -1.007825 + (mz - 34.968853)*(2/3)
    mz_M2pCl = 34.968853 + (mz - 34.968853)*(2/3)
    mz_M2m2HpNa = 20.97412 + (mz - 34.968853)*(2/3)
    mz_M2m2HpNapHCOOH = 66.9796 + (mz - 34.968853)*(2/3)
    mz_M2m2HpK = 36.948058 + (mz - 34.968853)*(2/3)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])

def Solo_M3m2HpNapHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 66.979600)/3
    mz_Cl = 34.968853 + (mz - 66.979600)/3
    mz_m2HpNa = 20.97412 + (mz - 66.979600)/3
    mz_m2HpNapHCOOH = 66.9796 + (mz - 66.979600)/3
    mz_m2HpK = 36.948058 + (mz - 66.979600)/3

    mz_M2mH = -1.007825 + (mz - 66.979600)*(2/3)
    mz_M2pCl = 34.968853 + (mz - 66.979600)*(2/3)
    mz_M2m2HpNa = 20.97412 + (mz - 66.979600)*(2/3)
    mz_M2m2HpNapHCOOH = 66.9796 + (mz - 66.979600)*(2/3)
    mz_M2m2HpK = 36.948058 + (mz - 66.979600)*(2/3)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])

def Solo_M3m2HpNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 20.97412)/3
    mz_Cl = 34.968853 + (mz - 20.97412)/3
    mz_m2HpNa = 20.97412 + (mz - 20.97412)/3
    mz_m2HpNapHCOOH = 66.9796 + (mz - 20.97412)/3
    mz_m2HpK = 36.948058 + (mz - 20.97412)/3

    mz_M2mH = -1.007825 + (mz - 20.97412)*(2/3)
    mz_M2pCl = 34.968853 + (mz - 20.97412)*(2/3)
    mz_M2m2HpNa = 20.97412 + (mz - 20.97412)*(2/3)
    mz_M2m2HpNapHCOOH = 66.9796 + (mz - 20.97412)*(2/3)
    mz_M2m2HpK = 36.948058 + (mz - 20.97412)*(2/3)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])

def Solo_M4mH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz + 1.007825)/4
    mz_Cl = 34.968853 + (mz + 1.007825)/4
    mz_m2HpNa = 20.97412 + (mz + 1.007825)/4
    mz_m2HpNapHCOOH = 66.9796 + (mz + 1.007825)/4
    mz_m2HpK = 36.948058 + (mz + 1.007825)/4

    mz_M2mH = -1.007825 + (mz + 1.007825)/2
    mz_M2pCl = 34.968853 + (mz + 1.007825)/2
    mz_M2m2HpNa = 20.97412 + (mz + 1.007825)/2
    mz_M2m2HpNapHCOOH = 66.9796 + (mz + 1.007825)/2
    mz_M2m2HpK = 36.948058 + (mz + 1.007825)/2

    mz_M3mH = -1.007825 + (mz + 1.007825)*(3/4)
    mz_M3pCl = 34.968853 + (mz + 1.007825)*(3/4)
    mz_M3m2HpNa = 20.97412 + (mz + 1.007825)*(3/4)
    mz_M3m2HpNapHCOOH = 66.9796 + (mz + 1.007825)*(3/4)
    mz_M3m2HpK = 36.948058 + (mz + 1.007825)*(3/4)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])

def Solo_M4pCl(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 34.968853)/4
    mz_Cl = 34.968853 + (mz - 34.968853)/4
    mz_m2HpNa = 20.97412 + (mz - 34.968853)/4
    mz_m2HpNapHCOOH = 66.9796 + (mz - 34.968853)/4
    mz_m2HpK = 36.948058 + (mz - 34.968853)/4

    mz_M2mH = -1.007825 + (mz - 34.968853)/2
    mz_M2pCl = 34.968853 + (mz - 34.968853)/2
    mz_M2m2HpNa = 20.97412 + (mz - 34.968853)/2
    mz_M2m2HpNapHCOOH = 66.9796 + (mz - 34.968853)/2
    mz_M2m2HpK = 36.948058 + (mz - 34.968853)/2

    mz_M3mH = -1.007825 + (mz - 34.968853)*(3/4)
    mz_M3pCl = 34.968853 + (mz - 34.968853)*(3/4)
    mz_M3m2HpNa = 20.97412 + (mz - 34.968853)*(3/4)
    mz_M3m2HpNapHCOOH = 66.9796 + (mz - 34.968853)*(3/4)
    mz_M3m2HpK = 36.948058 + (mz - 34.968853)*(3/4)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])

def Solo_M4m2HpNapHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 66.979600)/4
    mz_Cl = 34.968853 + (mz - 66.979600)/4
    mz_m2HpNa = 20.97412 + (mz - 66.979600)/4
    mz_m2HpNapHCOOH = 66.9796 + (mz - 66.979600)/4
    mz_m2HpK = 36.948058 + (mz - 66.979600)/4

    mz_M2mH = -1.007825 + (mz - 66.979600)/2
    mz_M2pCl = 34.968853 + (mz - 66.979600)/2
    mz_M2m2HpNa = 20.97412 + (mz - 66.979600)/2
    mz_M2m2HpNapHCOOH = 66.9796 + (mz - 66.979600)/2
    mz_M2m2HpK = 36.948058 + (mz - 66.979600)/2

    mz_M3mH = -1.007825 + (mz - 66.979600)*(3/4)
    mz_M3pCl = 34.968853 + (mz - 66.979600)*(3/4)
    mz_M3m2HpNa = 20.97412 + (mz - 66.979600)*(3/4)
    mz_M3m2HpNapHCOOH = 66.9796 + (mz - 66.979600)*(3/4)
    mz_M3m2HpK = 36.948058 + (mz - 66.979600)*(3/4)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])

def Solo_M4m2HpNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = -1.007825 + (mz - 20.97412)/4
    mz_Cl = 34.968853 + (mz - 20.97412)/4
    mz_m2HpNa = 20.97412 + (mz - 20.97412)/4
    mz_m2HpNapHCOOH = 66.9796 + (mz - 20.97412)/4
    mz_m2HpK = 36.948058 + (mz - 20.97412)/4

    mz_M2mH = -1.007825 + (mz - 20.97412)/2
    mz_M2pCl = 34.968853 + (mz - 20.97412)/2
    mz_M2m2HpNa = 20.97412 + (mz - 20.97412)/2
    mz_M2m2HpNapHCOOH = 66.9796 + (mz - 20.97412)/2
    mz_M2m2HpK = 36.948058 + (mz - 20.97412)/2

    mz_M3mH = -1.007825 + (mz - 20.97412)*(3/4)
    mz_M3pCl = 34.968853 + (mz - 20.97412)*(3/4)
    mz_M3m2HpNa = 20.97412 + (mz - 20.97412)*(3/4)
    mz_M3m2HpNapHCOOH = 66.9796 + (mz - 20.97412)*(3/4)
    mz_M3m2HpK = 36.948058 + (mz - 20.97412)*(3/4)
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])

def Solo_M2pH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 1.007825)/2
    mz_Na = 22.98977 + (mz - 1.007825)/2
    mz_K = 38.963708 + (mz - 1.007825)/2
    mz_HpCH3CN = 42.034374 + (mz - 1.007825)/2
    mz_HpCH3OH = 33.034040 + (mz - 1.007825)/2
    mz_NapCH3CN = 64.016319 + (mz - 1.007825)/2
    mz_NapCH3OH = 55.015985 + (mz - 1.007825)/2
    mz_KpCH3CN = 79.990257 + (mz - 1.007825)/2
    mz_KpCH3OH = 70.989923 + (mz - 1.007825)/2
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])

def Solo_M2pHpCH3CN(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 42.034374)/2
    mz_Na = 22.98977 + (mz - 42.034374)/2
    mz_K = 38.963708 + (mz - 42.034374)/2
    mz_HpCH3CN = 42.034374 + (mz - 42.034374)/2
    mz_HpCH3OH = 33.034040 + (mz - 42.034374)/2
    mz_NapCH3CN = 64.016319 + (mz - 42.034374)/2
    mz_NapCH3OH = 55.015985 + (mz - 42.034374)/2
    mz_KpCH3CN = 79.990257 + (mz - 42.034374)/2
    mz_KpCH3OH = 70.989923 + (mz - 42.034374)/2

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    
    return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])

def Solo_M2pHpCH3OH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 33.034040)/2
    mz_Na = 22.98977 + (mz - 33.034040)/2
    mz_K = 38.963708 + (mz - 33.034040)/2
    mz_HpCH3CN = 42.034374 + (mz - 33.034040)/2
    mz_HpCH3OH = 33.034040 + (mz - 33.034040)/2
    mz_NapCH3CN = 64.016319 + (mz - 33.034040)/2
    mz_NapCH3OH = 55.015985 + (mz - 33.034040)/2
    mz_KpCH3CN = 79.990257 + (mz - 33.034040)/2
    mz_KpCH3OH = 70.989923 + (mz - 33.034040)/2

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    
    return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])


def Solo_M2pHpHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 47.013304)/2
    mz_Na = 22.98977 + (mz - 47.013304)/2
    mz_K = 38.963708 + (mz - 47.013304)/2
    mz_HpCH3CN = 42.034374 + (mz - 47.0133042)/2
    mz_HpCH3OH = 33.034040 + (mz - 47.013304)/2
    mz_NapCH3CN = 64.016319 + (mz - 47.013304)/2
    mz_NapCH3OH = 55.015985 + (mz - 47.013304)/2
    mz_KpCH3CN = 79.990257 + (mz - 47.013304)/2
    mz_KpCH3OH = 70.989923 + (mz - 47.013304)/2

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])

def Solo_M2pNH4(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 18.034374)/2
    mz_NH4 = 18.034374 + (mz - 18.034374)/2
    mz_Na = 22.98977 + (mz - 18.034374)/2
    mz_K = 38.963708 + (mz - 18.034374)/2
    mz_HpCH3CN = 42.034374 + (mz - 18.034374)/2
    mz_HpCH3OH = 33.034040 + (mz - 18.034374)/2
    mz_NapCH3CN = 64.016319 + (mz - 18.034374)/2
    mz_NapCH3OH = 55.015985 + (mz - 18.034374)/2
    mz_KpCH3CN = 79.990257 + (mz - 18.034374)/2
    mz_KpCH3OH = 70.989923 + (mz - 18.034374)/2

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])

    
def Solo_M2pNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 22.98977)/2
    mz_Na = 22.98977 + (mz - 22.98977)/2
    mz_NapCH3CN = 64.016319 + (mz - 22.98977)/2
    mz_NapCH3OH = 55.015985 + (mz - 22.98977)/2
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])

def Solo_M2pNapCH3OH(prec_mass_error, ion_idx, spectrum_list) : 
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 55.015985)/2
    mz_Na = 22.98977 + (mz - 55.015985)/2
    mz_NapCH3CN = 64.016319 + (mz - 55.015985)/2
    mz_NapCH3OH = 55.015985 + (mz - 55.015985)/2
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])

def Solo_M2pNapCH3CN(prec_mass_error, ion_idx, spectrum_list) : 
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 64.016319)/2
    mz_Na = 22.98977 + (mz - 64.016319)/2
    mz_NapCH3CN = 64.016319 + (mz - 64.016319)/2
    mz_NapCH3OH = 55.015985 + (mz - 64.016319)/2
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])

def Solo_M2pK(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 38.963708)/2
    mz_Na = 22.98977 + (mz - 38.963708)/2
    mz_K = 38.963708 + (mz - 38.963708)/2

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na, valid_K])

def Solo_M1pHpCH3CN(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 42.034374
    mz_Na = 22.98977 + mz - 42.034374
    mz_HpCH3OH = 33.034040 + mz - 42.034374
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Na, valid_HpCH3OH])

def Solo_M1pHpCH3OH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 33.034040
    mz_Na = 22.98977 + mz - 33.034040
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Na])

def Solo_M1pHpHCOOH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 47.013304
    mz_Na = 22.98977 + mz - 47.013304
    mz_HpCH3OH = 33.034040 + mz - 47.013304
    mz_HpCH3CN = 42.034374 + mz - 47.013304
 
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0

    return sum([valid_H, valid_Na, valid_HpCH3OH, valid_HpCH3CN])

def Solo_M1pNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 22.989770

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()

    return valid_H

def Solo_M1pNapCH3CN(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = mz - 64.016319 + 1.007825
    mz_Na = mz - 64.016319 + 22.98977
    mz_NapCH3OH = mz - 64.016319 + 55.015985
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na, valid_NapCH3OH])

def Solo_M1pNapCH3OH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 55.015985
    mz_Na = 22.98977 + mz - 55.015985
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_Na])

def Solo_M1pNH4(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 18.034374
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
    return valid_H

def Solo_M1pNH4pCH3CN(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 59.060923
    mz_NH4 = 18.034374 + mz - 59.060923
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    return sum([valid_H, valid_NH4])

def Solo_M1pK(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 38.963708
    mz_NH4 =  18.034374 + mz - 38.963708
    mz_Na = 22.98977 + mz - 38.963708
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum()
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum()
    return sum([valid_H, valid_NH4, valid_Na])

def Solo_M1pKpCH3OH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + mz - 70.989923
    mz_NH4 =  18.034374 + mz - 70.989923
    mz_Na = 22.98977 + mz - 70.989923
    mz_K = 38.963708 + mz - 70.989923
    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum()
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum()
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum()
    return sum([valid_H, valid_NH4, valid_Na, valid_K])

def Solo_M3pH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 1.007825)/3
    mz_NH4 = 18.034374 + (mz - 1.007825)/3
    mz_Na = 22.98977 + (mz - 1.007825)/3
    mz_K = 38.963708 + (mz - 1.007825)/3
    mz_HpCH3CN = 42.034374 + (mz - 1.007825)/3
    mz_HpCH3OH = 33.034040 + (mz - 1.007825)/3
    mz_NapCH3CN = 64.016319 + (mz - 1.007825)/3
    mz_NapCH3OH = 55.015985 + (mz - 1.007825)/3
    mz_KpCH3CN = 79.990257 + (mz - 1.007825)/3
    mz_KpCH3OH = 70.989923 + (mz - 1.007825)/3
    mz_M2pH = 1.007825 + (mz - 1.007825)*(2/3)
    mz_M2pNH4 = 18.034374 + (mz - 18.034374)*(2/3)
    mz_M2pNa = 22.98977 + (mz - 1.007825)*(2/3)
    mz_M2pK = 38.963708 + (mz - 1.007825)*(2/3)
    mz_M2pHpCH3CN = 42.034374 + (mz - 1.007825)*(2/3)
    mz_M2pHpCH3OH = 33.034040 + (mz - 1.007825)*(2/3)
    mz_M2pNapCH3CN = 64.016319 + (mz - 1.007825)*(2/3)
    mz_M2pNapCH3OH = 55.015985 + (mz - 1.007825)*(2/3)
    mz_M2pKpCH3CN = 79.990257 + (mz - 1.007825)*(2/3)
    mz_M2pKpCH3OH = 70.989923 + (mz - 1.007825)*(2/3)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH])

def Solo_M3pNH4(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 18.034374)/3
    mz_NH4 = 18.034374 + (mz - 18.034374)/3
    mz_Na = 22.98977 + (mz - 18.034374)/3
    mz_K = 38.963708 + (mz - 18.034374)/3
    mz_HpCH3CN = 42.034374 + (mz - 18.034374)/3
    mz_HpCH3OH = 33.034040 + (mz - 18.034374)/3
    mz_NapCH3CN = 64.016319 + (mz - 18.034374)/3
    mz_NapCH3OH = 55.015985 + (mz - 18.034374)/3
    mz_KpCH3CN = 79.990257 + (mz - 18.034374)/3
    mz_KpCH3OH = 70.989923 + (mz - 18.034374)/3
    mz_M2pH = 1.007825 + (mz - 18.034374)*(2/3)
    mz_M2pNH4 = 18.034374 + (mz - 18.034374)*(2/3)
    mz_M2pNa = 22.98977 + (mz - 18.034374)*(2/3)
    mz_M2pK = 38.963708 + (mz - 18.034374)*(2/3)
    mz_M2pHpCH3CN = 42.034374 + (mz - 18.034374)*(2/3)
    mz_M2pHpCH3OH = 33.034040 + (mz - 18.034374)*(2/3)
    mz_M2pNapCH3CN = 64.016319 + (mz - 18.034374)*(2/3)
    mz_M2pNapCH3OH = 55.015985 + (mz - 18.034374)*(2/3)
    mz_M2pKpCH3CN = 79.990257 + (mz - 18.034374)*(2/3)
    mz_M2pKpCH3OH = 70.989923 + (mz - 18.034374)*(2/3)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH])

    
def Solo_M3pNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 22.989770)/3
    mz_NH4 = 18.034374 + (mz - 22.989770)/3
    mz_Na = 22.98977 + (mz - 22.989770)/3
    mz_K = 38.963708 + (mz - 22.989770)/3
    mz_HpCH3CN = 42.034374 + (mz - 22.989770)/3
    mz_HpCH3OH = 33.034040 + (mz - 22.989770)/3
    mz_NapCH3CN = 64.016319 + (mz - 22.989770)/3
    mz_NapCH3OH = 55.015985 + (mz - 22.989770)/3
    mz_KpCH3CN = 79.990257 + (mz - 22.989770)/3
    mz_KpCH3OH = 70.989923 + (mz - 22.989770)/3
    mz_M2pH = 1.007825 + (mz - 22.989770)*(2/3)
    mz_M2pNH4 = 18.034374 + (mz - 22.989770)*(2/3)
    mz_M2pNa = 22.98977 + (mz - 22.989770)*(2/3)
    mz_M2pK = 38.963708 + (mz - 22.989770)*(2/3)
    mz_M2pHpCH3CN = 42.034374 + (mz - 22.989770)*(2/3)
    mz_M2pHpCH3OH = 33.034040 + (mz - 22.989770)*(2/3)
    mz_M2pNapCH3CN = 64.016319 + (mz - 22.989770)*(2/3)
    mz_M2pNapCH3OH = 55.015985 + (mz - 22.989770)*(2/3)
    mz_M2pKpCH3CN = 79.990257 + (mz - 22.989770)*(2/3)
    mz_M2pKpCH3OH = 70.989923 + (mz - 22.989770)*(2/3)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH])

def Solo_M3pK(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 38.963708)/3
    mz_NH4 = 18.034374 + (mz - 38.963708)/3
    mz_Na = 22.98977 + (mz - 38.963708)/3
    mz_K = 38.963708 + (mz - 38.963708)/3
    mz_HpCH3CN = 42.034374 + (mz - 38.963708)/3
    mz_HpCH3OH = 33.034040 + (mz - 38.9637080)/3
    mz_NapCH3CN = 64.016319 + (mz - 38.9637080)/3
    mz_NapCH3OH = 55.015985 + (mz - 38.963708)/3
    mz_KpCH3CN = 79.990257 + (mz - 38.963708)/3
    mz_KpCH3OH = 70.989923 + (mz - 38.963708)/3
    mz_M2pH = 1.007825 + (mz - 38.963708)*(2/3)
    mz_M2pNH4 = 18.034374 + (mz - 38.963708)*(2/3)
    mz_M2pNa = 22.98977 + (mz - 38.963708)*(2/3)
    mz_M2pK = 38.963708 + (mz - 38.963708)*(2/3)
    mz_M2pHpCH3CN = 42.034374 + (mz - 38.963708)*(2/3)
    mz_M2pHpCH3OH = 33.034040 + (mz - 38.963708)*(2/3)
    mz_M2pNapCH3CN = 64.016319 + (mz - 38.963708)*(2/3)
    mz_M2pNapCH3OH = 55.015985 + (mz - 38.963708)*(2/3)
    mz_M2pKpCH3CN = 79.990257 + (mz - 38.963708)*(2/3)
    mz_M2pKpCH3OH = 70.989923 + (mz - 38.963708)*(2/3)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH])
    
def Solo_M4pK(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 38.963708)/4
    mz_NH4 = 18.034374 + (mz - 38.963708)/4
    mz_Na = 22.98977 + (mz - 38.963708)/4
    mz_K = 38.963708 + (mz - 38.963708)/4
    mz_HpCH3CN = 42.034374 + (mz - 38.963708)/4
    mz_HpCH3OH = 33.034040 + (mz - 38.963708)/4
    mz_NapCH3CN = 64.016319 + (mz - 38.963708)/4
    mz_NapCH3OH = 55.015985 + (mz - 38.963708)/4
    mz_KpCH3CN = 79.990257 + (mz - 38.963708)/4
    mz_KpCH3OH = 70.989923 + (mz - 38.963708)/4
    mz_M2pH = 1.007825 + (mz - 38.963708)/2
    mz_M2pNH4 = 18.034374 + (mz - 38.963708)/2
    mz_M2pNa = 22.98977 + (mz - 38.963708)/2
    mz_M2pK = 38.963708 + (mz - 38.963708)/2
    mz_M2pHpCH3CN = 42.034374 + (mz - 38.963708)/2
    mz_M2pHpCH3OH = 33.034040 + (mz - 38.963708)/2
    mz_M2pNapCH3CN = 64.016319 + (mz - 38.963708)/2
    mz_M2pNapCH3OH = 55.015985 + (mz - 38.963708)/2
    mz_M2pKpCH3CN = 79.990257 + (mz - 38.963708)/2
    mz_M2pKpCH3OH = 70.989923 + (mz - 38.963708)/2
    mz_M3pH = 1.007825 + (mz - 38.963708)*(3/4)
    mz_M3pNH4 = 18.034374 + (mz - 38.963708)*(3/4)
    mz_M3pNa = 22.98977 + (mz - 38.963708)*(3/4)
    mz_M3pK = 38.963708 + (mz - 38.963708)*(3/4)
    mz_M3pHpCH3CN = 42.034374 + (mz - 38.963708)*(3/4)
    mz_M3pHpCH3OH = 33.034040 + (mz - 38.963708)*(3/4)
    mz_M3pNapCH3CN = 64.016319 + (mz - 38.963708)*(3/4)
    mz_M3pNapCH3OH = 55.015985 + (mz - 38.963708)*(3/4)
    mz_M3pKpCH3CN = 79.990257 + (mz - 38.963708)*(3/4)
    mz_M3pKpCH3OH = 70.989923 + (mz - 38.963708)*(3/4)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                valid_M3pKpCH3CN, valid_M3pKpCH3OH])

def Solo_M4pH(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 1.007825)/4
    mz_NH4 = 18.034374 + (mz - 1.007825)/4
    mz_Na = 22.98977 + (mz - 1.007825)/4
    mz_K = 38.963708 + (mz - 1.007825)/4
    mz_HpCH3CN = 42.034374 + (mz - 1.007825)/4
    mz_HpCH3OH = 33.034040 + (mz - 1.007825)/4
    mz_NapCH3CN = 64.016319 + (mz - 1.007825)/4
    mz_NapCH3OH = 55.015985 + (mz - 1.007825)/4
    mz_KpCH3CN = 79.990257 + (mz - 1.007825)/4
    mz_KpCH3OH = 70.989923 + (mz - 1.007825)/4
    mz_M2pH = 1.007825 + (mz - 1.007825)/2
    mz_M2pNH4 = 18.034374 + (mz - 1.007825)/2
    mz_M2pNa = 22.98977 + (mz - 1.007825)/2
    mz_M2pK = 38.963708 + (mz - 1.007825)/2
    mz_M2pHpCH3CN = 42.034374 + (mz - 1.007825)/2
    mz_M2pHpCH3OH = 33.034040 + (mz - 1.007825)/2
    mz_M2pNapCH3CN = 64.016319 + (mz - 1.007825)/2
    mz_M2pNapCH3OH = 55.015985 + (mz - 1.007825)/2
    mz_M2pKpCH3CN = 79.990257 + (mz - 1.007825)/2
    mz_M2pKpCH3OH = 70.989923 + (mz - 1.007825)/2
    mz_M3pH = 1.007825 + (mz - 1.007825)*(3/4)
    mz_M3pNH4 = 18.034374 + (mz - 1.007825)*(3/4)
    mz_M3pNa = 22.98977 + (mz - 1.007825)*(3/4)
    mz_M3pK = 38.963708 + (mz - 1.007825)*(3/4)
    mz_M3pHpCH3CN = 42.034374 + (mz - 1.007825)*(3/4)
    mz_M3pHpCH3OH = 33.034040 + (mz - 1.007825)*(3/4)
    mz_M3pNapCH3CN = 64.016319 + (mz - 1.007825)*(3/4)
    mz_M3pNapCH3OH = 55.015985 + (mz - 1.007825)*(3/4)
    mz_M3pKpCH3CN = 79.990257 + (mz - 1.007825)*(3/4)
    mz_M3pKpCH3OH = 70.989923 + (mz - 1.007825)*(3/4)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                valid_M3pKpCH3CN, valid_M3pKpCH3OH])

def Solo_M4pNH4(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 18.034374)/4
    mz_NH4 = 18.034374 + (mz - 18.034374)/4
    mz_Na = 22.98977 + (mz - 18.034374)/4
    mz_K = 38.963708 + (mz - 18.034374)/4
    mz_HpCH3CN = 42.034374 + (mz - 18.034374)/4
    mz_HpCH3OH = 33.034040 + (mz - 18.034374)/4
    mz_NapCH3CN = 64.016319 + (mz - 18.034374)/4
    mz_NapCH3OH = 55.015985 + (mz - 18.034374)/4
    mz_KpCH3CN = 79.990257 + (mz - 18.034374)/4
    mz_KpCH3OH = 70.989923 + (mz - 18.034374)/4
    mz_M2pH = 1.007825 + (mz - 18.034374)/2
    mz_M2pNH4 = 18.034374 + (mz - 18.034374)/2
    mz_M2pNa = 22.98977 + (mz - 18.034374)/2
    mz_M2pK = 38.963708 + (mz - 18.034374)/2
    mz_M2pHpCH3CN = 42.034374 + (mz - 18.034374)/2
    mz_M2pHpCH3OH = 33.034040 + (mz - 18.034374)/2
    mz_M2pNapCH3CN = 64.016319 + (mz - 18.034374)/2
    mz_M2pNapCH3OH = 55.015985 + (mz - 18.034374)/2
    mz_M2pKpCH3CN = 79.990257 + (mz - 18.034374)/2
    mz_M2pKpCH3OH = 70.989923 + (mz - 18.034374)/2
    mz_M3pH = 1.007825 + (mz - 18.034374)*(3/4)
    mz_M3pNH4 = 18.034374 + (mz - 18.034374)*(3/4)
    mz_M3pNa = 22.98977 + (mz - 18.034374)*(3/4)
    mz_M3pK = 38.963708 + (mz - 18.034374)*(3/4)
    mz_M3pHpCH3CN = 42.034374 + (mz - 18.034374)*(3/4)
    mz_M3pHpCH3OH = 33.034040 + (mz - 18.034374)*(3/4)
    mz_M3pNapCH3CN = 64.016319 + (mz - 18.034374)*(3/4)
    mz_M3pNapCH3OH = 55.015985 + (mz - 18.034374)*(3/4)
    mz_M3pKpCH3CN = 79.990257 + (mz - 18.034374)*(3/4)
    mz_M3pKpCH3OH = 70.989923 + (mz - 18.034374)*(3/4)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0

    valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = True).sum() > 0
    valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = True).sum() > 0
    valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = True).sum() > 0
    valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = True).sum() > 0
    valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                valid_M3pKpCH3CN, valid_M3pKpCH3OH])

def Solo_M4pNa(prec_mass_error, ion_idx, spectrum_list):
    mz = spectrum_list[ion_idx].get('pepmass')[0]
    peaks = pd.Series(spectrum_list[ion_idx].peaks.mz)
    mz_H = 1.007825 + (mz - 22.98977)/4
    mz_NH4 = 18.034374 + (mz - 22.98977)/4
    mz_Na = 22.98977 + (mz - 22.98977)/4
    mz_K = 38.963708 + (mz - 22.98977)/4
    mz_HpCH3CN = 42.034374 + (mz - 22.98977)/4
    mz_HpCH3OH = 33.034040 + (mz - 22.98977)/4
    mz_NapCH3CN = 64.016319 + (mz - 22.98977)/4
    mz_NapCH3OH = 55.015985 + (mz - 22.98977)/4
    mz_KpCH3CN = 79.990257 + (mz - 22.98977)/4
    mz_KpCH3OH = 70.989923 + (mz - 22.98977)/4
    mz_M2pH = 1.007825 + (mz - 22.98977)/2
    mz_M2pNH4 = 18.034374 + (mz - 22.98977)/2
    mz_M2pNa = 22.98977 + (mz - 22.98977)/2
    mz_M2pK = 38.963708 + (mz - 22.98977)/2
    mz_M2pHpCH3CN = 42.034374 + (mz - 22.98977)/2
    mz_M2pHpCH3OH = 33.034040 + (mz - 22.98977)/2
    mz_M2pNapCH3CN = 64.016319 + (mz - 22.98977)/2
    mz_M2pNapCH3OH = 55.015985 + (mz - 22.98977)/2
    mz_M2pKpCH3CN = 79.990257 + (mz - 22.98977)/2
    mz_M2pKpCH3OH = 70.989923 + (mz - 22.98977)/2
    mz_M3pH = 1.007825 + (mz - 22.98977)*(3/4)
    mz_M3pNH4 = 18.034374 + (mz - 22.98977)*(3/4)
    mz_M3pNa = 22.98977 + (mz - 22.98977)*(3/4)
    mz_M3pK = 38.963708 + (mz - 22.98977)*(3/4)
    mz_M3pHpCH3CN = 42.034374 + (mz - 22.98977)*(3/4)
    mz_M3pHpCH3OH = 33.034040 + (mz - 22.98977)*(3/4)
    mz_M3pNapCH3CN = 64.016319 + (mz - 22.98977)*(3/4)
    mz_M3pNapCH3OH = 55.015985 + (mz - 22.98977)*(3/4)
    mz_M3pKpCH3CN = 79.990257 + (mz - 22.98977)*(3/4)
    mz_M3pKpCH3OH = 70.989923 + (mz - 22.98977)*(3/4)

    valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0

    valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        
    return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                valid_M3pKpCH3CN, valid_M3pKpCH3OH])
    
def Not_validable(prec_mass_error, ion_idx, spectrum_list):
    return 0

def Validator_choice(adduct, ion_mode):
    """Selects the Species_rule function appropriate to the provided adduct.
    """
    if ion_mode == "POS":
        if adduct == 'M1|p1H|pCH3OH' : return Solo_M1pHpCH3OH
        elif adduct == 'M1|p1H|pHCOOH' : return Solo_M1pHpHCOOH
        elif adduct == 'M2|p1Na|' : return Solo_M2pNa
        elif adduct == 'M1|p1NH4|' : return Solo_M1pNH4
        elif adduct == 'M1|p1NH4|pCH3CN' : return Solo_M1pNH4pCH3CN
        elif adduct == 'M1|p1Na|' : return Solo_M1pNa
        elif adduct == 'M1|p1H|' : return Not_validable
        elif adduct == 'M1|p1Na|pCH3CN' : return Solo_M1pNapCH3CN
        elif adduct == 'M1|p1H|pCH3CN' : return Solo_M1pHpCH3CN
        elif adduct == 'M2|p1H|' : return Solo_M2pH
        elif adduct == 'M2|p1H|pCH3CN' : return Solo_M2pHpCH3CN
        elif adduct == 'M2|p1NH4|' : return Solo_M2pNH4
        elif adduct == 'M2|p1K|' : return Solo_M2pK
        elif adduct == 'M1|p1K|' : return Solo_M1pK
        elif adduct == 'M1|p1K|pCH3OH' : return Solo_M1pKpCH3OH
        elif adduct == 'M2|p1Na|pCH3OH' : return Solo_M2pNapCH3OH
        elif adduct == 'M2|p1H|pHCOOH' : return Solo_M2pHpHCOOH
        elif adduct == 'M2|p1Na|pCH3CN' : return Solo_M2pNapCH3CN
        elif adduct == 'M2|p1H|pCH3OH' : return Solo_M2pHpCH3OH
        elif adduct == 'M3|p1H|' : return Solo_M3pH
        elif adduct == 'M3|p1Na|' : return Solo_M3pNa
        elif adduct == 'M3|p1NH4|' : return Solo_M3pNH4
        elif adduct == 'M4|p1K|' : return Solo_M4pK
        elif adduct == 'M4|p1H|' : return Solo_M4pH
        elif adduct == 'M4|p1NH4|' : return Solo_M4pNH4
        elif adduct == 'M4|p1Na|' : return Solo_M4pNa
        elif adduct == 'M3|p1K|' : return Solo_M3pK
    elif ion_mode == "NEG":
        if adduct == 'M1|m1H|' : return Not_validable
        elif adduct == 'M1|p1Cl|' : return Not_validable
        elif adduct == 'M1|m1H|pC4H11N' : return Solo_M1mHpC4H11N
        elif adduct == 'M1|m1H|pHCOOH' : return Solo_M1mHpHCOOH
        elif adduct == 'M1|m2Hp1Na|pHCOOH' : return Solo_M1m2HpNapHCOOH
        elif adduct == 'M1|m2Hp1Na|' : return Solo_M1m2HpNa
        elif adduct == 'M1|m2Hp1K|' : return Solo_M1m2HpK
        elif adduct == 'M2|m1H|pC4H11N' : return Solo_M2mHpC4H11N
        elif adduct == 'M2|m1H|pHCOOH' : return Solo_M2mHpHCOOH
        elif adduct == 'M2|m1H|' : return Solo_M2mH
        elif adduct == 'M2|p1Cl|' : return Solo_M2pCl
        elif adduct == 'M2|m2Hp1Na|pHCOOH' : return Solo_M2m2HpNapHCOOH
        elif adduct == 'M2|m2Hp1Na|' : return Solo_M2m2HpNa
        elif adduct == 'M2|m2Hp1K|' : return Solo_M2m2HpK
        elif adduct == 'M3|m1H|' : return Solo_M3mH
        elif adduct == 'M3|p1Cl|' : return Solo_M3pCl
        elif adduct == 'M3|m2Hp1Na|pHCOOH' : return Solo_M3m2HpNapHCOOH
        elif adduct == 'M3|m2Hp1Na|' : return Solo_M3m2HpNa
        elif adduct == 'M4|m1H|' : return Solo_M4mH
        elif adduct == 'M4|p1Cl|' : return Solo_M4pCl
        elif adduct == 'M4|m2Hp1Na|pHCOOH' : return Solo_M4m2HpNapHCOOH
        elif adduct == 'M4|m2Hp1Na|' : return Solo_M4m2HpNa