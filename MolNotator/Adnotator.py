def Adnotator(params : dict, ion_mode : str):
    
    def Spectrum_processing(s):
        s = default_filters(s)
        return s
    
    def Solo_M1mHpC4H11N(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 72.081324
        mz_Cl = 34.968853 + mz - 72.081324
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl])
    
    def Solo_M1mHpHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 44.997654
        mz_Cl = 34.968853 + mz - 44.997654
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl])
    
    def Solo_M1m2HpNapHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 66.979600
        mz_Cl = 34.968853 + mz - 66.979600
        mz_m2HpNa = 20.97412 + mz - 66.979600
        
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])


    def Solo_M1m2HpNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 66.979600
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H])
    
    def Solo_M1m2HpK(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 36.948058
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H])
    
    def Solo_M2mHpC4H11N(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 72.081324)/2
        mz_Cl = 34.968853 + (mz - 72.081324)/2
        mz_m2HpNa = 20.97412 + (mz - 72.081324)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])
    
    def Solo_M2mHpHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 44.997654)/2
        mz_Cl = 34.968853 + (mz - 44.997654)/2
        mz_m2HpNa = 20.97412 + (mz - 44.997654)/2
        mz_mHpHCOOH = 44.997654 + (mz - 44.997654)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
        valid_mHpHCOOH = peaks.between(mz_mHpHCOOH - prec_mass_error, mz_mHpHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_mHpHCOOH])
    
    def Solo_M2mH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz + 1.007825)/2
        mz_Cl = 34.968853 + (mz + 1.007825)/2
        mz_m2HpNa = 20.97412 + (mz + 1.007825)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])
    
    def Solo_M2pCl(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 34.968853)/2
        mz_Cl = 34.968853 + (mz - 34.968853)/2
        mz_m2HpNa = 20.97412 + (mz - 34.968853)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])
    
    def Solo_M2m2HpNapHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 66.979600)/2
        mz_Cl = 34.968853 + (mz - 66.979600)/2
        mz_m2HpNa = 20.97412 + (mz - 66.979600)/2
        mz_m2HpNapHCOOH = 66.9796 + (mz - 66.979600)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH])
    
    def Solo_M2m2HpNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 20.97412)/2
        mz_Cl = 34.968853 + (mz - 20.97412)/2
        mz_m2HpNa = 20.97412 + (mz - 20.97412)/2
        mz_m2HpNapHCOOH = 66.9796 + (mz - 20.97412)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = "both").sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH])
    
    def Solo_M2m2HpK(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M3mH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M3pCl(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M3m2HpNapHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M3m2HpNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4mH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4pCl(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4m2HpNapHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4m2HpNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M2pH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M2pHpCH3CN(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M2pHpCH3OH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    
    def Solo_M2pHpHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M2pNH4(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
        
    def Solo_M2pNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 22.98977)/2
        mz_Na = 22.98977 + (mz - 22.98977)/2
        mz_NapCH3CN = 64.016319 + (mz - 22.98977)/2
        mz_NapCH3OH = 55.015985 + (mz - 22.98977)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])
    
    def Solo_M2pNapCH3OH(ion_idx) : 
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 55.015985)/2
        mz_Na = 22.98977 + (mz - 55.015985)/2
        mz_NapCH3CN = 64.016319 + (mz - 55.015985)/2
        mz_NapCH3OH = 55.015985 + (mz - 55.015985)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])
    
    def Solo_M2pNapCH3CN(ion_idx) : 
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 64.016319)/2
        mz_Na = 22.98977 + (mz - 64.016319)/2
        mz_NapCH3CN = 64.016319 + (mz - 64.016319)/2
        mz_NapCH3OH = 55.015985 + (mz - 64.016319)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])
    
    def Solo_M2pK(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 38.963708)/2
        mz_Na = 22.98977 + (mz - 38.963708)/2
        mz_K = 38.963708 + (mz - 38.963708)/2
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_Na, valid_K])
    
    def Solo_M1pHpCH3CN(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 42.034374
        mz_Na = 22.98977 + mz - 42.034374
        mz_HpCH3OH = 33.034040 + mz - 42.034374
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Na, valid_HpCH3OH])
    
    def Solo_M1pHpCH3OH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 33.034040
        mz_Na = 22.98977 + mz - 33.034040
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Na])
    
    def Solo_M1pHpHCOOH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 47.013304
        mz_Na = 22.98977 + mz - 47.013304
        mz_HpCH3OH = 33.034040 + mz - 47.013304
        mz_HpCH3CN = 42.034374 + mz - 47.013304
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = "both").sum() > 0
    
        return sum([valid_H, valid_Na, valid_HpCH3OH, valid_HpCH3CN])
    
    def Solo_M1pNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 22.989770
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
    
        return valid_H
    
    def Solo_M1pNapCH3CN(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = mz - 64.016319 + 1.007825
        mz_Na = mz - 64.016319 + 22.98977
        mz_NapCH3OH = mz - 64.016319 + 55.015985
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3OH])
    
    def Solo_M1pNapCH3OH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 55.015985
        mz_Na = 22.98977 + mz - 55.015985
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_Na])
    
    def Solo_M1pNH4(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 18.034374
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
        return valid_H
    
    def Solo_M1pNH4pCH3CN(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 59.060923
        mz_NH4 = 18.034374 + mz - 59.060923
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum() > 0
        return sum([valid_H, valid_NH4])
    
    def Solo_M1pK(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 38.963708
        mz_NH4 =  18.034374 + mz - 38.963708
        mz_Na = 22.98977 + mz - 38.963708
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum()
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum()
        return sum([valid_H, valid_NH4, valid_Na])
    
    def Solo_M1pKpCH3OH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 70.989923
        mz_NH4 =  18.034374 + mz - 70.989923
        mz_Na = 22.98977 + mz - 70.989923
        mz_K = 38.963708 + mz - 70.989923
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = "both").sum()
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = "both").sum()
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = "both").sum()
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = "both").sum()
        return sum([valid_H, valid_NH4, valid_Na, valid_K])
    
    def Solo_M3pH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M3pNH4(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
        
    def Solo_M3pNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M3pK(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
        
    def Solo_M4pK(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4pH(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4pNH4(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
    
    def Solo_M4pNa(ion_idx):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
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
        
    def Not_validable(adduct):
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
    def Get_feature_ids(full_mgf_file):
        from matchms.importing import load_from_mgf
        full_mgf_file = list(load_from_mgf(full_mgf_file))
        feature_ids = []
        for spectrum in full_mgf_file:
            feature_ids.append(int(spectrum.get("feature_id")))
        feature_ids = pd.Series(index = feature_ids, data = list(range(len(feature_ids))))
        return feature_ids
    
    def Cross_sample_tables(full_mgf_file):
        """Create dataframe to report results using the complete original MGF file.
        Each dataframe will contain data for ions (indexes) for each given sample 
        (columns). 
        """
        from matchms.importing import load_from_mgf
        full_mgf_file = list(load_from_mgf(full_mgf_file))
        feature_ids = []
        for spectrum in full_mgf_file:
            feature_ids.append(int(spectrum.get("feature_id")))
        cross_annotations = pd.DataFrame(index = feature_ids, dtype = str)
        cross_points = pd.DataFrame(index = feature_ids, dtype = float)
        cross_courts = pd.DataFrame(index = feature_ids, dtype = int)
        cross_houses = pd.DataFrame(index = feature_ids, dtype = int)
        cross_rules = pd.DataFrame(index = feature_ids, dtype = int)
        cross_neutrals = pd.DataFrame(index = feature_ids, dtype = int)
        feature_ids = pd.Series(index = feature_ids, data = list(range(len(feature_ids))))
        return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals, feature_ids
    
    def Rt_slicer(rt, rt_error, ion_id, input_table) :
        """Returns a dataframe containing all ions coeluted with the target ion (ion_id)
        given its rt (retention time) and the rt_error supplied by the user
        """
        rt_low = rt - rt_error
        rt_high = rt + rt_error
        sliced_table = input_table[input_table['rt'].between(rt_low,
                                  rt_high, inclusive = "both")].copy()
        return sliced_table.drop(ion_id)
    
    def Charge_filter(max_charge, input_table):
        """Eliminates ions which charge above max_charge. For now, Adnotator can
        only deal with single charge ions, max_charge is set at 1 by default
        """
        return input_table[abs(input_table['charge']) <= max_charge]
    
    def Neutral_to_adduct(mol_mass, adduct_mass, mol_count, ion_charge) :
        """Calculates the m/z value of an ion species, given its supposed molecular
        mass (neutral) and other parameters from the adducts table.
        """
        return round((mol_count*mol_mass + adduct_mass)/abs(ion_charge), 4)
    
    def Neutral_mass_calculator(input_mz, adduct_mass, mol_count, ion_charge) :
        """Calculates the molecular mass of a molecule given an ion species, its 
        m/z value and other parameters from the adducts table.
        """
        return round(((input_mz*abs(ion_charge)) - adduct_mass)/mol_count, 4)
    
    def Neutral_table(ion1_mz, adduct_table_primary):
        """Computes all possible molecular masses (hypothetical neutrals) for an ion
        (ion 1) given the different ion species available in the adducts table.
        """
        neutral_table = adduct_table_primary.copy()
        neutral_table['neutral_mass'] = [0.0]*len(neutral_table)
        for i in neutral_table.index :
            neutral_table.loc[i, 'neutral_mass'] = Neutral_mass_calculator(ion1_mz, 
                             neutral_table['Adduct_mass'][i], 
                             neutral_table['Mol_multiplier'][i], 
                             neutral_table['Charge'][i])
        return neutral_table
    
    def Ion_hypotheses_table(neutral_table):
        """Computes all possible ion m/z values (ion 2) given all possible molecular masses
        computed before with the Neutral_table function (hypothetical neutrals), 
        for a single ion (ion 1).
        """
        ion1_adduct = list()
        ion2_adduct = list()
        ion2_complexities = list()
        ion2_mz = list()
        for i in neutral_table.index :
            for j in neutral_table.index :
                if i == j: continue
                ion1_adduct.append(i) 
                ion2_adduct.append(j)
                ion2_complexities.append(neutral_table['Complexity'][j])
                ion2_mz.append(Neutral_to_adduct(neutral_table['neutral_mass'][i],
                                                  neutral_table['Adduct_mass'][j],
                                                  neutral_table['Mol_multiplier'][j],
                                                  neutral_table['Charge'][j]))
        ion_hypotheses_table = list(zip(ion1_adduct, ion2_adduct, ion2_complexities, ion2_mz))
        return pd.DataFrame(ion_hypotheses_table, columns = ['Ion1_adduct', 'Ion2_adduct', 'Ion2_complexity', 'Ion2_mz'])
    
    def Point_counter(ion_hypotheses_table, mass_error):
        """Checks among the coeluted ions (in the coelution table) if any of them
        matches an "ion 2" from the ion_hypotheses_table. Eliminates all ion 2 
        hypotheses that could not be matched to any coeluted ion.
        """
        hits = list()
        ion_indexes = list()
        for i in ion_hypotheses_table.index :
            ion2_mz = ion_hypotheses_table['Ion2_mz'][i]
            ion2_mz_low = ion2_mz - mass_error
            ion2_mz_high = ion2_mz + mass_error
            hit_count = coelution_table['mz'].between(ion2_mz_low, ion2_mz_high, 
                                      inclusive = "both").sum()
            if hit_count > 0 :
                hit_ids = list(coelution_table.index[coelution_table['mz'].between(ion2_mz_low, 
                                            ion2_mz_high, inclusive = "both")])
            else:
                hit_ids = []
            hits.append(hit_count)
            ion_indexes.append(hit_ids)
        ion_hypotheses_table['hit_count'] = hits
        ion_hypotheses_table['hit_indexes'] = ion_indexes
        return ion_hypotheses_table[ion_hypotheses_table['hit_count'] > 0]
            
    
    def Fragnotator_points(ion1_idx, ion_hypotheses_table, edge_table, status):
        """Awards points for each hypothesis in which the ions 1 and 2 are already
        linked by a fragmentation bond by fragnotator.
        """
        ion_hypotheses_table = ion_hypotheses_table.copy()
        if status != "fragment":
            temp_edges = edge_table[edge_table['node_1'] == ion1_idx]
            if len(temp_edges) == 0 : return ion_hypotheses_table
            frag_bond = list()
            for ion2 in ion_hypotheses_table['hit_indexes']:
                if ion2 in list(temp_edges['node_2']):
                    frag_bond.append(1)
                else :
                    frag_bond.append(0)
            ion_hypotheses_table['fragnotator_points'] = frag_bond
        else :
            temp_edges = edge_table[edge_table['node_2'] == ion1_idx]
            if len(temp_edges) == 0 : return ion_hypotheses_table
            frag_bond = list()
            for ion2 in ion_hypotheses_table['hit_indexes']:
                if ion2 in list(temp_edges['node_1']):
                    frag_bond.append(1)
                else :
                    frag_bond.append(0)
            ion_hypotheses_table['fragnotator_points'] = frag_bond
        return ion_hypotheses_table
        
    def Complex_points(neutral_table, ion_hypotheses_table):
        """For a hypothesis bearing an ion 2 with a complex form (neutral complexes,
        i.e. the third componend of the adduct code : Mol|ion|complex), points are
        awared if the uncomplexed form is found among the other hypotheses.
        """
        ion_hypotheses_table['Complex_points'] = [0]*len(ion_hypotheses_table)
        for i in ion_hypotheses_table.index :
            ion1_idx = ion_hypotheses_table['Ion1_adduct'][i]
            
            # Is the level 2 ion hypothesis a solvent complex?
            # Get the level 2 ion hypothesis code
            ion2_code = neutral_table.loc[ion_hypotheses_table['Ion2_adduct'][i], 'Adduct_code']
            # Check if there is indeed a solvent complex value present
            ion2_code_split = ion2_code.split('|')
            if ion2_code_split[2] != "" :
                # Get the mass of the uncomplexed form :
                ion2_code_split[2] = ""
                ion2_uncomplexed = '|'.join(ion2_code_split)
                ion2_uncomplexed_idx = neutral_table.index[neutral_table['Adduct_code']==ion2_uncomplexed][0]
                
                # Check if the uncomplexed form was detected using the level_2 table filtered
                bool_1 = ion_hypotheses_table['Ion1_adduct'] == ion1_idx
                bool_2 = ion_hypotheses_table['Ion2_adduct'] == ion2_uncomplexed_idx
                
                # Reward with points level 2 ion hypotheses for which uncomplexed 
                # ions could be detected
                if (bool_1 & bool_2).sum() > 0 : 
                    ion_hypotheses_table.loc[i, 'Complex_points'] += 1
        return ion_hypotheses_table
    
    def Get_adduct_code(ion_hypotheses_table, neutral_table):
        """Adds the adduct code for each hypothesis' ion 1.
        """
        ion_hypotheses_table['adduct_code'] = ['']*len(ion_hypotheses_table)
        for i in ion_hypotheses_table.index :
            ion_hypotheses_table.loc[i, 'adduct_code'] = neutral_table['Adduct_code'][ion_hypotheses_table['Ion1_adduct'][i]]
        return ion_hypotheses_table
    
    def Weighted_points(ion_hypotheses_table):
        """Calculated the points accumulated by each hypothesis, considering cosine 
        similarity, rule points, complex points, fragnotator points and the ion 2
        complexity.
        """
        ion_hypotheses_table['weighted_points'] = [0.0]*len(ion_hypotheses_table)
        for i in ion_hypotheses_table.index :
            ion_hypotheses_table.loc[i, 'weighted_points'] = ((ion_hypotheses_table['hit_count'][i] + 
                                    ion_hypotheses_table['cosine_score'][i]+
                                    ion_hypotheses_table['rule_points'][i]+
                                    ion_hypotheses_table['Complex_points'][i]+
                                    ion_hypotheses_table['fragnotator_points'][i])/
                                    (ion_hypotheses_table['Ion2_complexity'][i]+
                                     adduct_table_primary.loc[ion_hypotheses_table.loc[i, 'Ion1_adduct'], "Complexity"]))
        return ion_hypotheses_table
    
    def Merged_adducts_table(ion_hypotheses_table, ion_1_index):
        """Merges the ion_hypotheses_table for the current ion (i) being processed
        to all previously processed ions from the sample.
        """
        if len(ion_hypotheses_table) == 0 : return(ion_hypotheses_table)
        unique_adducts = list(ion_hypotheses_table['adduct_code'].unique())
        merged_adducts_table = pd.DataFrame(index = range(len(unique_adducts)))
        merged_adducts_table['ion_id'] = [ion_1_index]*len(unique_adducts)
        merged_adducts_table['adduct'] = unique_adducts
        merged_adducts_table['relation_codes'] = ['']*len(merged_adducts_table)
        merged_adducts_table['hit_count'] = [0]*len(merged_adducts_table)
        merged_adducts_table['hit_indexes'] = ['']*len(merged_adducts_table)
        merged_adducts_table['weighted_points'] = [0.0]*len(merged_adducts_table)
        for i in merged_adducts_table.index :
            adduct = merged_adducts_table.loc[i, 'adduct']
            tmp_table_1 = ion_hypotheses_table[ion_hypotheses_table['adduct_code'] == adduct]
            ion1_adduct = tmp_table_1['Ion1_adduct'].astype(str)
            ion2_adduct = tmp_table_1['Ion2_adduct'].astype(str)
            merged_adducts_table.loc[i, 'relation_codes'] = '|'.join(ion1_adduct+":"+ion2_adduct)
            merged_adducts_table.loc[i, 'hit_count'] = tmp_table_1['hit_count'].sum()
            merged_adducts_table.loc[i, 'hit_indexes'] = '|'.join(tmp_table_1['hit_indexes'].astype(str))
            merged_adducts_table.loc[i, 'weighted_points'] = tmp_table_1['weighted_points'].sum()
        return merged_adducts_table
    
    def Cohort_Table(neutral_table, full_merged_table, node_table):
        """Produce a cohort_table, with the ion_indexes as rows (ion IDs) and the
        ionisation hypotheses indexes as columns (taken from the full_merged_table).
        Each ion (row) will have 0, 1 or several ionisation hypotheses (columns). 
        and each of these columns will affect different ions based on the number 
        of ion_2 ionisation hypotheses that were found. Related ionisation hypotheses
        (ion ID + adduct formula) that are not mutually exclusive (i.e. different 
        ionisations for a single ion) are named "Cohorts".
        """
        ion_indexes = list(node_table.index)
        cohort_table = pd.DataFrame(data = None, index = ion_indexes,
                                    columns = full_merged_table.index)
        print('Producing cohort table...')
        for i in tqdm(full_merged_table.index) :
            current_ion = int(full_merged_table['ion_id'][i])
            split_relations_code = full_merged_table.loc[i,'relation_codes'].split('|')
            split_hit_ids = full_merged_table.loc[i,'hit_indexes'].split('|')
            ion1_adduct = int(split_relations_code[0].split(':')[0])
            ion1_adduct = neutral_table.loc[ion1_adduct, 'Adduct_code']
            cohort_table.loc[current_ion, i] = ion1_adduct
            for j in range(len(split_relations_code)) :
                split_relations_code[j] = split_relations_code[j].split(':')[1]
                split_relations_code[j] = neutral_table.loc[int(split_relations_code[j]),
                                    'Adduct_code']
            for j in range(len(split_hit_ids)) :
                split_ids = split_hit_ids[j]
                cohort_table.loc[int(split_ids), i] = split_relations_code[j]
        return cohort_table
    
    def Get_Court_Table(cohort_table) :
        """Produces a court_table from a cohort_table, containing all Courts produced
        from the different cohorts and the associated ions. A Court is a collection
        of related ion IDs and cohorts. Two cohorts, each with their own ions, 
        sharing at least one ion, will be considered related.
        """
        def Get_Court(cohort_list, cohort_table):
            """Produces a court from one or several cohorts. A Court is a collection
            of related ion IDs and cohorts. Two cohorts, each with their own ions,
            sharing at least one ion, will be considered related.
            """
            def Get_cohort_ions(cohorts, cohort_table) :
                """Returns every ion ID associated to the submitted cohorts list
                """
                if type(cohorts) != list:
                    cohorts = [cohorts]
                ion_list = list()
                for i in cohorts :
                    ion_list += list(cohort_table[i].dropna().index)
                ion_list = list(set(ion_list))
                ion_list.sort()
                return ion_list
            
            def Get_ions_cohorts(ion_list, cohort_table):
                """Returns every cohorts associated with the submitted ions list
                """
                if type(ion_list) != list:
                    ion_list = [ion_list]
                cohort_list = list()
                for i in ion_list :
                    cohort_list += list(cohort_table.loc[i].dropna().index)
                cohort_list = list(set(cohort_list))
                cohort_list.sort()
                return cohort_list
    
            nb_ions = 0
            nb_cohorts = 1
            court_expanding = True
            while court_expanding :
                ion_list = Get_cohort_ions(cohort_list, cohort_table)
                cohort_list = Get_ions_cohorts(ion_list, cohort_table)
                if nb_ions + nb_cohorts == len(ion_list) + len(cohort_list) :
                    court_expanding = False
                else :
                    nb_ions = len(ion_list)
                    nb_cohorts = len(cohort_list)
            return ion_list, cohort_list
        
        cohort_table_to_court = cohort_table.copy() # surrogate cohort_table to be emptied
        court_id = 0 # court_id represe
        court_table = pd.DataFrame(columns = ['Cohort_List', 'Ion_List'])
        while (len(cohort_table_to_court.index) > 0 and len(cohort_table_to_court.columns) > 0):
            # i is the cohort number to be analysed, to determine the entire court
            # it is involved with
            i = cohort_table_to_court.columns[0]
            # Get all related ions and cohorts and store them in the court_table
            ion_list, cohort_list = Get_Court(i, cohort_table_to_court)
            court_table.loc[court_id] = [cohort_list, ion_list]
            
            # remove all involved ions and cohorts from the cohort_to_court_table
            cohort_table_to_court.drop(ion_list, axis = 0, inplace = True)
            cohort_table_to_court.drop(cohort_list, axis = 1, inplace = True)
            court_id += 1
        # Add to the court_table all singleton ions (no cohorts associated)
        for i in cohort_table_to_court.index :
            court_table.loc[court_id] = [[], [i]]
            court_id += 1
        if len(court_table.index) > 0:
            court_table['Cohort_count'] = [0]*len(court_table)
            court_table['Ion_count'] = [0]*len(court_table)
            for i in court_table.index :
                court_table.loc[i, 'Cohort_count'] = len(court_table.loc[i, 'Cohort_List'])
                court_table.loc[i, 'Ion_count'] = len(court_table.loc[i, 'Ion_List'])
        court_table.sort_values("Cohort_count", ascending = False, inplace = True)
        court_table = court_table[court_table['Cohort_count'] > 0]
        court_table.reset_index(drop = True, inplace = True)
        return court_table
    
    def Cosine_validation(ion_1_idx, cosine_threshold, mgf_file):
        """ For each ionisation hypothesis (ion 1 and 2 pairs), reports the
        cosine score from the cosine_table. In case of multiple hits in the
        "hit_indexes" columns, the one with the best cosine score is kept (chosen
        at random if both have the same score : either both completely unrelated
        or the same spectrum). The duplicate hits are transferred to a duplicate_df
        for annotation in the final network, they will take no part in calculation.
        Points are then awared according to the cosine score:
        1+cosine score if the score is above threshold, 0 if it is below.
        """
        cosine_list = list()
        matched_peaks_list = list()
        product_list = list()
        different_group = list()

        ion_hypotheses_table.reset_index(drop = True, inplace = True)
        duplicate_df = pd.DataFrame(columns = ["ion_1_idx", "ion_1_adduct",
                                               "ion_2_idx", "ion_2_adduct",
                                               "selected_ion", "cosine_score",
                                               "matched_peaks", "different_group"])
        for i in ion_hypotheses_table.index:
            group_1 = adduct_table_primary.loc[ion_hypotheses_table.loc[i, "Ion1_adduct"], "Group"]
            group_2 = adduct_table_primary.loc[ion_hypotheses_table.loc[i, "Ion2_adduct"], "Group"]
            if group_1 == group_2: different_group.append(False)
            else: different_group.append(True)
            if ion_hypotheses_table.loc[i, "hit_count"] == 1 :
                ion_2_idx = ion_hypotheses_table.loc[i,"hit_indexes"][0]
                score, n_matches = modified_cosine.pair(mgf_file[ion_1_idx],
                                                        mgf_file[ion_2_idx])
                prod = score * n_matches
            else:
                selected_hit = pd.DataFrame(columns = ['cos', 'peaks', 'prod'])
                for hit in ion_hypotheses_table.loc[i,"hit_indexes"]:
                    score, n_matches = modified_cosine.pair(mgf_file[ion_1_idx],
                                                            mgf_file[hit])
                    selected_hit.loc[hit] = [score, n_matches, score * n_matches]

                selected_hit.sort_values('prod', ascending = False, inplace =True)
                new_hit = selected_hit.index[0]
                score = selected_hit['cos'].iloc[0]
                n_matches = int(selected_hit['peaks'].iloc[0])
                prod = selected_hit['prod'].iloc[0]
                selected_hit.drop(new_hit, inplace = True)
                ion_hypotheses_table.loc[i, "hit_indexes"] = [[new_hit]]
                for j in selected_hit.index:
                    tmp_idx = len(duplicate_df)
                    duplicate_df.loc[tmp_idx] = [ion_1_idx,
                                        ion_hypotheses_table.loc[i, "Ion1_adduct"],
                                     j, ion_hypotheses_table.loc[i, "Ion2_adduct"],
                                     new_hit, selected_hit.loc[j, 'cos'], selected_hit.loc[j, 'peaks'],
                                     different_group[-1]]
            if score < cosine_threshold: 
                score = 0.0
                prod = 0.0
            cosine_list.append(score)
            matched_peaks_list.append(n_matches)
            product_list.append(prod)
        ion_hypotheses_table['different_group'] = different_group
        ion_hypotheses_table['cosine_score'] = cosine_list
        ion_hypotheses_table['matched_peaks'] = matched_peaks_list
        ion_hypotheses_table['product'] = product_list
        
        ion_hypotheses_table['hit_indexes'] = ion_hypotheses_table['hit_indexes'].str[0]
        tmp_bool = ion_hypotheses_table["different_group"] + ion_hypotheses_table["cosine_score"] >= cosine_threshold
        duplicate_bool = duplicate_df["different_group"] + duplicate_df["cosine_score"] >= cosine_threshold
        return ion_hypotheses_table[tmp_bool], duplicate_df[duplicate_bool]
    
    def Species_rules(ion_1_idx, ion_hypotheses_table):
        """Awards points to hypotheses for which the annotations for ions 1 or 2 
        could be validated by the ion species rules (checks if the annotation can 
        be confirmed).
        """
        ion_hypotheses_table['rule_points'] = [0]*len(ion_hypotheses_table)
        for hypothesis in ion_hypotheses_table.index:
            ion_1_adduct = ion_hypotheses_table.loc[hypothesis, "Ion1_adduct"]
            ion_1_adduct = adduct_table_primary.loc[ion_1_adduct, "Adduct_code"]
            Species_rule = Validator_choice(ion_1_adduct, ion_mode)
            ion_hypotheses_table.loc[hypothesis, "rule_points"] += Species_rule(ion_1_idx)
            
            ion_2_adduct = ion_hypotheses_table.loc[hypothesis, "Ion2_adduct"]
            ion_2_adduct = adduct_table_primary.loc[ion_2_adduct, "Adduct_code"]
            ion_2_idx = ion_hypotheses_table.loc[hypothesis, "hit_indexes"]
            Species_rule = Validator_choice(ion_2_adduct, ion_mode)
            ion_hypotheses_table.loc[hypothesis, "rule_points"] += Species_rule(ion_2_idx)
        return ion_hypotheses_table
    
    def Cohort_ions(cohort):
        """Gets ions for a given supercohort (from the supercohorts table)
        """
        return list(supercohorts_table[cohort].dropna().index)
    def Incompatible_cohorts(cohort, cohort_pool):
        """Gets supercohorts incompatible with the given supercohort.
        """
        cohort_ions = Cohort_ions(cohort)
        incompatibles = []
        for ion in cohort_ions:
            annotation = supercohorts_table.loc[ion, cohort]
            incompatible = supercohorts_table.loc[ion].dropna() != annotation
            incompatibles += list(incompatible.index[incompatible])
        return list(set(incompatibles))
    
    def House_selection(court_table):
        """Finds the houses for each given court, i.e. supercohort combinations
        in a court that are compatible with each other. Once all combinations are 
        explored, the most likely house is each cohort is selected based on the points
        cumulated by their cohorts. If there is a draw between two or more houses, they
        are all kept. Final selection will happen when cross examining samples.
        """
        selected_houses = []
        house_count = []
        house_points = []
        unresolved_courts = list(court_table.index[court_table['Cohort_count'] > 1])
        resolved_courts = list(court_table.index[court_table['Cohort_count'] <= 1])
    
        for court_id in unresolved_courts:
            houses = []
            court_cohorts = court_table.loc[court_id, "Cohort_List"].copy()
            ion_list = list()
            for cohort in court_cohorts :
                ion_list += supercohorts_table[cohort].dropna().index.tolist()
            ion_list = list(set(ion_list))
            ion_list.sort()
            cohorts_table = pd.DataFrame(data = None, index = court_cohorts, columns = ion_list)
            for i in cohorts_table.index:
                for j in cohorts_table.columns:
                    cohorts_table.loc[i, j] = supercohorts_table.loc[j, i]
            # Get contested ions and conflict table
            conflict_table = list()
            for i in cohorts_table.columns:
                conflicted = cohorts_table[i].dropna().index.tolist()
                if len(conflicted) > 1 :
                    conflict_table.append((i, conflicted))
            conflict_table = pd.DataFrame(conflict_table, columns = ['conflicted_ion', 'conflicted_cohorts'])
            migration = list()
            is_mirror = list()
            for i in conflict_table.index:
                conflict_ion = conflict_table.loc[i, "conflicted_ion"]
                cohort_selection = list()
                for j in conflict_table.loc[i, "conflicted_cohorts"]:
                    cohort_ions = cohorts_table.loc[j].dropna().index.tolist()
                    cohort_ions.remove(conflict_ion)
                    score_list = list()
                    for k in cohort_ions:
                        score, n_matches = modified_cosine.pair(mgf_file[conflict_ion],
                                                                mgf_file[k])
                        if score <= cosine_threshold : score = 0
                        score_list.append((score, n_matches, score*n_matches))
                    score_list = pd.DataFrame(score_list, columns = ['score', 'matches', 'product'])
                    score_list.sort_values('product', ascending = False, inplace = True)
                    cohort_selection.append((j,
                                             score_list['score'].iloc[0],
                                             score_list['matches'].iloc[0],
                                             score_list['product'].iloc[0]))
                cohort_selection = pd.DataFrame(cohort_selection, columns = ['cohort', 'cos', 'matches', 'product'])
                cohort_selection.sort_values('product', ascending = False, inplace = True)
                migration_blocks = list()
                for prod in cohort_selection['product'].unique():
                    migration_blocks.append(list(cohort_selection['cohort'][cohort_selection['product'] == prod]))
                migration.append(migration_blocks)
        
                if len(migration_blocks) == 1 :
                    is_mirror.append(True)
                else:
                    is_mirror.append(False)
            conflict_table['migration'] = migration
            conflict_table['is_mirror'] = is_mirror
            
            migrating = True
            while migrating :
                migration_table = pd.DataFrame(columns = ['ion_list', 'ion_count'])
                for i in cohorts_table.index:
                    migration_table.loc[i] = [cohorts_table.loc[i].dropna().index.tolist(), len(cohorts_table.loc[i].dropna())]
                for i in conflict_table.index:
                    if conflict_table.loc[i, "is_mirror"] : continue
                    conflict_ion = conflict_table.loc[i, "conflicted_ion"]
                    looser_cohorts = list(flatten(conflict_table.loc[i, "migration"][1:]))
                    for j in looser_cohorts:
                        migration_table.loc[j, "ion_list"].remove(conflict_ion)
                        migration_table.loc[j, "ion_count"] = len(migration_table.loc[j, "ion_list"])
                dead_cohorts = migration_table.index[migration_table['ion_count'] <= 1]
                cohorts_table.drop(dead_cohorts, inplace = True)
                if len(dead_cohorts) == 0 : 
                    migrating = False
                    continue
                for i in conflict_table.index:
                    conflict_table.at[i, "conflicted_cohorts"] = list(set(conflict_table.loc[i, "conflicted_cohorts"]) - set(dead_cohorts))
                    tmp_migration = list()
                    for j in range(len(conflict_table.loc[i, "migration"])):
                        tmp_cohorts = conflict_table.loc[i, "migration"][j]
                        tmp_cohorts = list(set(tmp_cohorts) - set(dead_cohorts))
                        if len(tmp_cohorts) > 0 :
                            tmp_migration.append(tmp_cohorts)
                    conflict_table.at[i, "migration"] = tmp_migration
            # Check mirrors or no-choice migrations again:
            for i in conflict_table.index:
                if len(conflict_table.loc[i, "migration"]) == 1 :
                    conflict_table.loc[i, "is_mirror"] = True
            
            # Do migrations :
            for i in conflict_table.index:
                if conflict_table.loc[i, "is_mirror"] : continue
                loosing_cohorts = list(flatten(conflict_table.loc[i, "migration"][1:]))
                conflict_ion = conflict_table.loc[i, "conflicted_ion"]
                for j in loosing_cohorts:
                    cohorts_table.loc[j,conflict_ion] = None
            
            # Check for dead cohorts due to duplicated annotations:
            dead_cohorts = list()
            for i in cohorts_table.index:
                if len(cohorts_table.loc[i].dropna().unique()) <= 1 : dead_cohorts.append(i)
            cohorts_table.drop(dead_cohorts, inplace = True)
            
            # Correct migrations in supercohorts table
            for i in cohorts_table.index:
                cohort_ions = cohorts_table.loc[i].dropna().index.tolist()
                other_ions = list(set(supercohorts_table[i].dropna().index.tolist()) - set(cohort_ions))
                for j in other_ions:
                    supercohorts_table.loc[j, i] = None
            
            # Check for compatibilities, save full compatible cohorts, proceed to house selection for the others.
            inert_cohorts = list()
            for i in cohorts_table.index:
                cohort_ions = cohorts_table.loc[i].dropna().index.tolist()
                full_compatible = True
                for ion in cohort_ions:
                    if len(cohorts_table[ion].dropna().unique()) > 1 : 
                        full_compatible = False
                        break
                if full_compatible:
                    inert_cohorts.append(i)
            cohorts_table.drop(inert_cohorts, inplace = True)

            court_cohorts = cohorts_table.index.tolist()
            while len(court_cohorts) > 0 :
                seed_cohort = court_cohorts[0]
                court_cohorts.remove(seed_cohort)
                # Initialise house table:
                next_list, house_list, compatibles_list, left_list = Initialise_house_table(seed_cohort, court_cohorts)
                if len(house_list) == 0 : houses.append([seed_cohort])
                while len(house_list) > 0 :
                    next_list, house_list, compatibles_list, left_list, houses = Clean_house_table(next_list, house_list, compatibles_list, left_list, houses)
                    next_list, house_list, compatibles_list, left_list = Expand_house_table(next_list, house_list, compatibles_list, left_list)
            
            # Get points for these cohorts
            house_table = pd.DataFrame(index = range(len(houses)))
            house_table['house'] = houses
            tmp_house_points = []
            tmp_house_ions = []
            print('Assigning points to each house...')
            for i in tqdm(house_table.index):
                tmp_points = 0.0
                tmp_ion_count = []
                for supercohort in house_table.loc[i, "house"]:
                    #supercohort = house_table.loc[i, "house"][0]
                    tmp_ion_count += list(supercohorts_table[supercohort].dropna().index)
                    cohorts = transition_table.loc[supercohort, "cohorts"]
                    tmp_points_list = []
                    for cohort in cohorts :
                        tmp_points_list.append(full_merged_table.loc[cohort, "weighted_points"])
                    tmp_points += max(tmp_points_list)
                tmp_ion_count = list(set(tmp_ion_count))
                tmp_house_ions.append(len(tmp_ion_count))
                tmp_house_points.append(tmp_points)
            house_table['points'] = tmp_house_points
            house_table['ion_count'] = tmp_house_ions
            house_table.sort_values(by = ['points', 'ion_count'], ascending = False, inplace = True)
            
            if len(house_table) > 0 :
                for i in inert_cohorts:
                    tmp_points = []
                    tmp_ion_count = len(supercohorts_table[i].dropna())
                    cohort_list = transition_table.loc[i, "cohorts"]
                    for cohort in cohort_list:
                        tmp_points.append(full_merged_table.loc[cohort, "weighted_points"])
                    tmp_points = max(tmp_points)
                    house_table['points'] += tmp_points
                    for j in house_table.index:
                        house_table.at[j, "house"] += inert_cohorts
                    house_table['ion_count'] += tmp_ion_count
                for i in house_table.index:
                    house_table.at[i, "house"] = [house_table.loc[i, "house"]]
            else :
                tmp_points = []
                tmp_ion_count = 0
                for i in inert_cohorts:
                    cohort_list = transition_table.loc[i, "cohorts"]
                    tmp_ion_count += len(cohort_list)
                    for cohort in cohort_list :
                        tmp_points.append(full_merged_table.loc[cohort, "weighted_points"])
                tmp_points = max(tmp_points)
                house_table = pd.DataFrame(columns = ['house', 'points', 'ion_count'])
                house_table.loc[0] = [[inert_cohorts], tmp_points, tmp_ion_count]
            selected_houses.append(house_table['house'].iloc[0])
            house_count.append(len(house_table['house'].iloc[0]))
            house_points.append(house_table['points'].max())
            
        for court_id in resolved_courts:
            selected_houses.append([court_table.loc[court_id, "Cohort_List"]])
            house_count.append(1)
            supercohort = court_table.loc[court_id, "Cohort_List"][0]
            cohorts = transition_table.loc[supercohort, "cohorts"]
            tmp_points = []
            for cohort in cohorts:
                tmp_points.append(full_merged_table.loc[cohort, "weighted_points"])
            house_points.append(max(tmp_points))
        court_table["selected_houses"] = selected_houses
        court_table["house_count"] = house_count
        court_table["house_points"] = house_points
        return court_table
    
    def Initialise_house_table(seed_cohort, court_cohorts):
        """Starts house building by combining compatible cohorts, using an initial
        cohort as a seed.
        """
        incompatibles = Incompatible_cohorts(seed_cohort, court_cohorts)
        compatibles = list(set(court_cohorts) - set(incompatibles))
        next_list = compatibles
        house_list = [[seed_cohort]]*len(next_list)
        compatibles_list = [None]*len(next_list)
        left_list = [0]*len(next_list)
        added = []
        for i in range(len(next_list)):
            added.append(next_list[i])
            tmp_incompatibles = Incompatible_cohorts(next_list[i], court_cohorts)
            tmp_incompatibles += added
            compatibles_list[i] = list(set(compatibles) - set(tmp_incompatibles))
            left_list[i] = len(compatibles_list[i])
        return next_list, house_list, compatibles_list, left_list
        
        # Clean house table and export results
    def Clean_house_table(next_list, house_list, compatibles_list, left_list, houses):
        """Removes houses that finished expanding (no compatible cohorts left) from
        the house table, and exports them to the houses list.
        """
        export_idx = [i for i in range(len(left_list)) if left_list[i] == 0]
        export_idx.sort(reverse = True)
        for i in export_idx :
            new_house = house_list[i] + [next_list[i]]
            new_house = list(set(new_house))
            new_house.sort()
            houses.append(new_house) 
            next_list.pop(i)
            house_list.pop(i)
            compatibles_list.pop(i)
            left_list.pop(i)
        return next_list, house_list, compatibles_list, left_list, houses
    
    def Expand_house_table(next_list, house_list, compatibles_list, left_list):
        """Expands the houses in the house table by appending to them one of each
        possible compatible cohorts.
        """
        drop_list = []
        print("Expanding house table...")
        for i in tqdm(range(len(next_list))):
            drop_list.append(i)
            new_house = house_list[i] + [next_list[i]]
            added = []
            compatibles = compatibles_list[i]
            for j in compatibles:
                #j = compatibles[0]
                added.append(j)
                incompatibles = Incompatible_cohorts(j, compatibles)
                incompatibles += added
                next_list.append(j)
                house_list.append(new_house)
                compatibles_list.append(list(set(compatibles) - set(incompatibles)))
                left_list.append(len(list(set(compatibles) - set(incompatibles))))
        drop_list.sort(reverse = True)
        for i in drop_list:
            next_list.pop(i)
            house_list.pop(i)
            compatibles_list.pop(i)
            left_list.pop(i)
        return next_list, house_list, compatibles_list, left_list
            
    def Supercohorts_table(cohort_table):
        """Merges compatible cohorts into a single supercohort to reduce combinations
        and speed up calculation times. Returns the supercohorts_table and the
        transitions_table which keeps track of the cohorts forming each supercohort.
        """
        # Super cohorts:
        transition_table = pd.DataFrame(columns = ["supercohort", "cohorts"])
        supercohorts_table = pd.DataFrame(index = cohort_table.index)
        cohort_pool = pd.Series(cohort_table.columns, index = cohort_table.columns)
        while len(cohort_pool) > 0 :
            i = cohort_pool.iloc[0]
            supercohort = cohort_table[i].dropna()
            mergeable_cohorts = [i]
            all_ions = list(supercohort.index)
            added_ions = all_ions
            len_0 = 0
            while len(all_ions) != len_0:
                len_0 = len(all_ions)
                for j in added_ions:
                    annotation = supercohort[j]
                    other_cohorts = cohort_table.loc[j].dropna()
                    other_cohorts = other_cohorts[list(set(other_cohorts.index) - set(mergeable_cohorts))]
                    other_cohorts = other_cohorts.index[other_cohorts == annotation]
                    mergeable_cohorts += list(other_cohorts)
                    added_ions = []
                    for k in other_cohorts :
                        new_ions = list(set(cohort_table[k].dropna().index) - set(all_ions))
                        added_ions += new_ions
                        for l in new_ions:
                            supercohort[l] = cohort_table.loc[l, k]
                    all_ions += added_ions
    
            if len(mergeable_cohorts) == 0 :
                supercohorts_table[len(supercohorts_table.columns)] = cohort_table[i]
                transition_table.loc[len(transition_table)] = [len(transition_table), [i]]
                cohort_pool.drop(i, inplace = True)
                continue
            mergeable_cohorts = list(set(mergeable_cohorts))
            mergeable_cohorts.sort()
    
            supercohorts_table[len(supercohorts_table.columns)] = [None]*len(supercohorts_table)
            for j in supercohort.index:
                supercohorts_table.loc[j, len(supercohorts_table.columns) -1] = supercohort[j]
            transition_table.loc[len(transition_table)] = [len(transition_table), mergeable_cohorts]
            cohort_pool.drop(mergeable_cohorts, inplace = True)
        return supercohorts_table, transition_table
    
    def Cross_sample_report(court_table):
        """Reports the data for the sample being processed into the cross tables 
        before exporting them as csv files.
        """
        cross_annotations[sample_base_name] = [None]*len(cross_annotations)
        cross_points[sample_base_name] = [None]*len(cross_points)
        cross_courts[sample_base_name] = [None]*len(cross_courts)
        cross_houses[sample_base_name] = [None]*len(cross_courts)
        cross_rules[sample_base_name] = [None]*len(cross_courts)
        cross_neutrals[sample_base_name] = [None]*len(cross_courts)
        print("Reporting results...")
        for i in tqdm(court_table.index):
            for house in range(len(court_table.loc[i, "selected_houses"])):
                cohorts = court_table.loc[i, "selected_houses"][house]
                for cohort in cohorts:
                    annotations = supercohorts_table[cohort].dropna()
                    for ion in annotations.index:
                        neutral_mass = adduct_table_merged["Adduct_code"] == annotations[ion]
                        neutral_mass = neutral_mass.index[neutral_mass][0]
                        neutral_mass = Neutral_mass_calculator(node_table.loc[ion, "mz"],
                                                               adduct_table_merged.loc[neutral_mass, "Adduct_mass"],
                                                               adduct_table_merged.loc[neutral_mass, "Mol_multiplier"],
                                                               adduct_table_merged.loc[neutral_mass, "Charge"])
                        Species_rule = Validator_choice(annotations[ion], ion_mode)
                        rule_validation = Species_rule(ion)
                        ion_id = int(node_table.loc[ion, "feature_id"])
                        if cross_annotations.loc[ion_id, sample_base_name] == None :
                            cross_annotations.loc[ion_id, sample_base_name] = annotations[ion]
                            cross_houses.loc[ion_id, sample_base_name] = str(house)
                            cross_rules.loc[ion_id, sample_base_name] = str(rule_validation)
                            cross_neutrals.loc[ion_id, sample_base_name] = str(neutral_mass)
                        else:
                            cross_annotations.loc[ion_id, sample_base_name] += "&" + annotations[ion]
                            cross_houses.loc[ion_id, sample_base_name] += "&" + str(house)
                            cross_rules.loc[ion_id, sample_base_name] += "&" + str(rule_validation)
                            cross_neutrals.loc[ion_id, sample_base_name] += "&" + str(neutral_mass)
                        cross_points.loc[ion_id, sample_base_name] = court_table.loc[i, "house_points"]
                        cross_courts.loc[ion_id, sample_base_name] = i
    
        # Add duplicate annotations :
        for i in duplicate_df.index:
            ion_2_idx = duplicate_df.loc[i, "ion_2_idx"]
            ion_2_id = int(node_table.loc[ion_2_idx, "feature_id"])
            if type(cross_annotations.loc[ion_2_id, sample_base_name]) == str : continue
            ion_2_adduct = duplicate_df.loc[i, "ion_2_adduct"]
            ion_2_adduct = adduct_table_primary.loc[ion_2_adduct, "Adduct_code"]
            selected_idx = duplicate_df.loc[i, "selected_ion"]
            selected_id = int(node_table.loc[selected_idx, "feature_id"])
            selected_adduct = cross_annotations.loc[selected_id, sample_base_name]
            if ion_2_adduct == selected_adduct:
                cross_annotations.loc[ion_2_id, sample_base_name] = ion_2_adduct
                cross_houses.loc[ion_2_id, sample_base_name] = cross_houses.loc[selected_id, sample_base_name]
                cross_points.loc[ion_2_id, sample_base_name] = cross_points.loc[selected_id, sample_base_name]
                cross_courts.loc[ion_2_id, sample_base_name] = cross_courts.loc[selected_id, sample_base_name]
                cross_neutrals.loc[ion_2_id, sample_base_name] = cross_neutrals.loc[selected_id, sample_base_name]
                cross_rules.loc[ion_2_id, sample_base_name] = cross_rules.loc[selected_id, sample_base_name]
        return
    
    def Secondary_adducts(mgf_file):
        """After processing a sample and annotating ions with the adducts in the 
        primary adducts table, the remaining ions are annotated with adducts from
        the secondary adducts table. This is done on the basis of the already
        computed molecules and consumes less resources. 
        """
        annotated = list(cross_neutrals[sample_base_name].dropna().index)
        unnannotated = list(set(node_table['feature_id'].astype(int)) - set(annotated))
        unnannotated = [node_table.index[node_table["feature_id"] == i][0] for i in unnannotated]
        unnannotated_table = node_table.loc[unnannotated]
        for ion_id in annotated:
            neutrals = cross_neutrals.loc[ion_id, sample_base_name].split('&')
            annotations = cross_annotations.loc[ion_id, sample_base_name].split('&')
            ion_idx = node_table.index[node_table["feature_id"] == ion_id][0]
            neutral_rt = node_table.loc[ion_idx, "rt"]
            coelution_table = unnannotated_table[unnannotated_table['rt'].between(neutral_rt - rt_error,  neutral_rt + rt_error, inclusive = "both")]
            houses = cross_houses.loc[ion_id, sample_base_name].split('&')
            for k in range(len(houses)): 
                house = houses[k]
                neutral = float(neutrals[int(k)])
                main_adduct = annotations[int(k)]
                main_group = adduct_table_primary["Group"][adduct_table_primary["Adduct_code"] == main_adduct].iloc[0]
                for i in adduct_table_secondary.index:
                    tmp_adduct = adduct_table_secondary.loc[i, "Adduct_code"]
                    tmp_mz= Neutral_to_adduct(neutral, adduct_table_secondary.loc[i, "Adduct_mass"],
                                              adduct_table_secondary.loc[i, "Mol_multiplier"],
                                              adduct_table_secondary.loc[i, "Charge"])
                    tmp_group = adduct_table_secondary.loc[i, "Group"]
                    tmp_table = coelution_table['mz'].between(tmp_mz - mass_error, tmp_mz + mass_error, inclusive = "both")
                    if sum(tmp_table) != 0:
                        #print(ion_id)
                        tmp_table = coelution_table[tmp_table]
                        for j in tmp_table.index:
                            valid = False
                            Species_rule = Validator_choice(tmp_adduct, ion_mode)
                            if main_group == tmp_group:
                                tmp_cos, n_matches = modified_cosine.pair(mgf_file[ion_idx],
                                                                        mgf_file[j])
                                if tmp_cos >= cosine_threshold : valid = True
                            elif Species_rule(j) != 0:
                                valid = True
                            if valid :
                                j_id = int(tmp_table.loc[j, "feature_id"])
                                cross_courts.loc[j_id, sample_base_name] = cross_courts.loc[ion_id, sample_base_name]
                                cross_points.loc[j_id, sample_base_name] = cross_points.loc[ion_id, sample_base_name]
                                
                                if cross_annotations.loc[j_id, sample_base_name] == None:
                                    cross_annotations.loc[j_id, sample_base_name] = tmp_adduct
                                    cross_houses.loc[j_id, sample_base_name] = house
                                    cross_neutrals.loc[j_id, sample_base_name] = str(neutral)
                                    cross_rules.loc[j_id, sample_base_name] = str(Species_rule(j))
                                elif tmp_adduct not in cross_annotations.loc[j_id, sample_base_name].split('&'):
                                    cross_annotations.loc[j_id, sample_base_name] += "&" + tmp_adduct
                                    cross_houses.loc[j_id, sample_base_name] += "&" + house
                                    cross_neutrals.loc[j_id, sample_base_name] += "&" + str(neutral)
                                    cross_rules.loc[j_id, sample_base_name] += "&" + str(Species_rule(j))
    ###############################################################################
    ########################## PART II FUNCTIONS ##################################
    ###############################################################################
    def Get_sample_houses(cross_houses, sample):
        sample_houses = cross_houses[sample].dropna().str.split('&')
        for i in sample_houses.index:
            sample_houses[i] = list(map(int, sample_houses[i]))
        return sample_houses
    
    def Get_sample_courts(court_table, sample):
        return court_table[sample].dropna()
    
    def Get_house_ions(sample_houses, house):
        return set([k for k in sample_houses.index if house in sample_houses[k]])
    
    def Get_court_ions(sample_courts, court):
        return set(sample_courts.index[sample_courts == court])
    
    def Get_cohort_ions(court_ions, house_ions):
        return list(court_ions.intersection(house_ions))
    
    
    
    def Resolver_level_2(ion_id, sample):
        sample_houses = Get_sample_houses(cross_houses, sample)
        sample_courts = Get_sample_courts(cross_courts, sample)
        conflict_table = pd.DataFrame()
        conflict_table['annotation'] = cross_annotations.loc[ion_id, sample].split('&')
        conflict_table['house_IDs'] = sample_houses[ion_id]
        conflict_table['neutrals'] = cross_neutrals.loc[ion_id, sample].split('&')
        conflict_table['rule_points'] = cross_rules.loc[ion_id, sample].split('&')
        conflict_table['conflicting_ions'] = [[]]*len(conflict_table)
        conflict_table['points'] = [0.0]*len(conflict_table)
        court = sample_courts[ion_id]
        court_ions = Get_court_ions(sample_courts, court)
        conflicting_ions = []
        for i in conflict_table.index:
            house = conflict_table.loc[i, "house_IDs"]
            house_ions = Get_house_ions(sample_houses, house)
            conflicting_ions.append(Get_cohort_ions(court_ions, house_ions))
        for i in conflict_table.index:
            house_ions = set(conflicting_ions[i])
            other_ions = conflicting_ions.copy()
            other_ions.remove(other_ions[i])
            other_ions = set(flatten(other_ions))
            conflict_table.at[i, "conflicting_ions"] = list(house_ions - other_ions)
        for i in conflict_table.index:
            for j in conflict_table.loc[i, "conflicting_ions"]:
                j_annotation = cross_annotations.loc[j, sample]
                if "&" in j_annotation : 
                    #continue #Try continue?
                    sys.exit(f"ERROR: conflict inside conflict for ion pair {ion_id}-{j}, sample {sample}")
                conflict_samples = cross_annotations.loc[j].dropna() == j_annotation
                conflict_samples = conflict_samples.index[conflict_samples]
                conflict_table.loc[i, "points"] += cross_points.loc[j, conflict_samples].max()
            
        conflict_table.sort_values("points", ascending = False, inplace = True)
        corrected_samples = cross_annotations.loc[ion_id] == cross_annotations.loc[ion_id, sample]
        corrected_samples = corrected_samples.index[corrected_samples]
        for i in corrected_samples:
            cross_annotations.loc[ion_id, i] = conflict_table["annotation"].iloc[0]
            cross_houses.loc[ion_id, i] = str(conflict_table["house_IDs"].iloc[0])
            cross_neutrals.loc[ion_id, i] = conflict_table["neutrals"].iloc[0]
            cross_rules.loc[ion_id, i] = conflict_table["rule_points"].iloc[0]
        conflict_table.drop(conflict_table.index[0], inplace = True)
        for i in conflict_table.index:
            for j in conflict_table.loc[i, "conflicting_ions"]:
                j_annotation = cross_annotations.loc[j, sample]
                corrected_samples = cross_annotations.loc[j] == j_annotation
                corrected_samples = corrected_samples.index[corrected_samples]
                for k in corrected_samples:
                    cross_annotations.loc[j, k] = None
                    cross_houses.loc[j, k] = None
                    cross_neutrals.loc[j, k] = None
                    cross_rules.loc[j, k] = None
                    cross_courts.loc[j, k] = None
                    cross_points.loc[j, k] = None
        return
    
    def Resolver_level_1(ion_id, sample):
        conflict_table = pd.DataFrame()
        conflict_table['annotation'] = cross_annotations.loc[ion_id, sample].split('&')
        conflict_table['house'] = cross_houses.loc[ion_id, sample].split('&')
        conflict_table['neutral'] = cross_neutrals.loc[ion_id, sample].split('&')
        conflict_table['rule_points'] = cross_rules.loc[ion_id, sample].split('&')
        conflict_table['points'] = [0.0]*len(conflict_table)
        for i in conflict_table.index:        
            tmp_samples = cross_annotations.loc[ion_id] == conflict_table.loc[i, "annotation"] 
            tmp_samples = tmp_samples.index[tmp_samples]
            conflict_table.loc[i, "points"] = cross_points.loc[ion_id, tmp_samples].max()
        selection_table = conflict_table[conflict_table['points'] == conflict_table['points'].max()]
        if len(selection_table) != 1 :
            Resolver_level_2(ion_id, sample)
            return
        else:
            conflict_table.drop(selection_table.index, inplace = True)
            correct_samples = cross_annotations.loc[ion_id] == cross_annotations.loc[ion_id, sample]
            correct_samples = correct_samples.index[correct_samples]
            for tmp_sample in correct_samples:
                cross_annotations.loc[ion_id, tmp_sample] = selection_table['annotation'].iloc[0]
                cross_houses.loc[ion_id, tmp_sample] = selection_table['house'].iloc[0]
                cross_neutrals.loc[ion_id, tmp_sample] = selection_table['neutral'].iloc[0]
                cross_rules.loc[ion_id, tmp_sample] = selection_table['rule_points'].iloc[0]
            for i in conflict_table.index:
                correct_samples = cross_annotations.loc[ion_id] == conflict_table.loc[i, "annotation"]
                correct_samples = correct_samples.index[correct_samples]
                for tmp_sample in correct_samples :
                    cross_annotations.loc[ion_id, tmp_sample] = None
                    cross_houses.loc[ion_id, tmp_sample] = None
                    cross_neutrals.loc[ion_id, tmp_sample] = None
                    cross_rules.loc[ion_id, tmp_sample] = None
                    cross_courts.loc[ion_id, tmp_sample] = None
                    cross_points.loc[ion_id, tmp_sample] = None
            return
    
    def Get_cohort_ions_2(ion, sample):
        sample_courts = Get_sample_courts(cross_courts, sample)
        sample_houses = Get_sample_houses(cross_houses, sample)
        court = cross_courts.loc[ion, sample]
        if court == None : return []
        house = int(cross_houses.loc[ion, sample])
        court_ions = Get_court_ions(sample_courts, court)
        house_ions = Get_house_ions(sample_houses, house)
        return Get_cohort_ions(court_ions, house_ions)
    
    def Merged_tables(input_files, feature_ids):
        """Runs through all samples in input_files (all user-submitted samples) and
        merges the fragnotator outputted tables into a single node and edge table,
        with the "status" updated for each ion in the node table.
        """
        merged_node_table = pd.DataFrame(columns = ['mz', 'rt', 'TIC', 'charge', 'mgf_index'])
        merged_edge_table = pd.DataFrame(columns = ['node_1', 'node_2', 'matched_peaks',
                                                    'total_peaks', 'matching_score',
                                                    'rt_gap', 'mz_gap', 'status',
                                                    'Fragnotation', 'All_annotations'])
        print('Merging node and edge tables from all samples:')
        for x in tqdm(input_files.index) :
            file_basename = input_files.at[x, 'base_name']    
            edge_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_edges.csv", 
                                     index_col = 'Index')
            node_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_nodes.csv", 
                                    index_col = 'Index')
            node_table['feature_id'] = node_table['feature_id'].astype(int)
            #Convert edge_table nodes to feature IDs
            for i in edge_table.index:
                edge_table.loc[i, 'node_1'] = node_table.loc[edge_table.loc[i, 'node_1'], "feature_id"]
                edge_table.loc[i, 'node_2'] = node_table.loc[edge_table.loc[i, 'node_2'], "feature_id"]
            # Retrieve mgf indexes from the feature IDs:
            tmp_mgf_idx = []
            for i in node_table.index:
                tmp_mgf_idx.append(feature_ids[node_table.loc[i, "feature_id"]])
            node_table.drop('status', axis = 1, inplace = True)
            node_table['mgf_index'] = tmp_mgf_idx
            node_table.set_index('feature_id', drop = True, inplace = True)
            new_ions = list(set(node_table.index) - set(merged_node_table.index))
            merged_node_table = merged_node_table.append(node_table.loc[new_ions], ignore_index = False)
            edge_table = edge_table[edge_table['status'] != "singleton"]
            edge_ids = edge_table['node_1'].astype(str) + "|" + edge_table['node_2'].astype(str)
            edge_table.set_index(edge_ids, inplace = True)
            new_edges = list(set(edge_table.index) - set(merged_edge_table.index))
            merged_edge_table = merged_edge_table.append(edge_table.loc[new_edges], ignore_index = False)
        merged_node_table.sort_index(inplace = True)
        merged_edge_table.sort_values("node_1", inplace = True)
        merged_edge_table.reset_index(drop = True, inplace = True)
        merged_node_table['status'] = ['singleton']*len(merged_node_table)
        precursors = list(set(merged_edge_table['node_1']) - set(merged_edge_table['node_2']))
        fragments = list(set(merged_edge_table['node_2']) - set(precursors))
        print('Updating precursor node status:')
        for i in tqdm(precursors):
            merged_node_table.loc[i, "status"] = "precursor"
        print('Updating fragment node status:')
        for i in tqdm(fragments):
            merged_node_table.loc[i, "status"] = "fragment"
        merged_edge_table = merged_edge_table.replace({np.nan: None})
        return merged_node_table, merged_edge_table
    
    def Annotation_resolver():
        """Resolves the yet sample-specific unresolved ions by cross examining 
        between all samples
        """
        print("Cross-examining unresolved annotations...")
        for sample in tqdm(cross_annotations.columns):
            tmp_annotations = cross_annotations[sample].dropna()
            tmp_annotations = tmp_annotations[tmp_annotations.str.contains('&')]
            if len(tmp_annotations) == 0 : continue
            for ion_id in tmp_annotations.index:
                Resolver_level_1(ion_id, sample)
        
        # Flush out annotations supported only by a single ion (after the above filtering)
        print("Flushing dead annotations...")
        for sample in tqdm(cross_annotations.columns):
            courts = cross_courts[sample].dropna()
            for court in courts.unique():
                court_samples = courts.index[courts==court]
                houses = cross_houses.loc[court_samples, sample]
                for house in houses.unique():
                    if (houses == house).sum() <= 1:
                        purge_index = houses[houses == house].index
                        for i in purge_index:
                            cross_annotations.loc[i, sample] = None
                            cross_courts.loc[i, sample] = None
                            cross_houses.loc[i, sample] = None
                            cross_neutrals.loc[i, sample] = None
                            cross_points.loc[i, sample] = None
                            cross_rules.loc[i, sample] = None
        return
    
    
    def Cross_court_table():
        """Produce a cross-sample court table to select the most likely houses
        valid across samples. compatible
        """
        print('Court processing - isolating singleton ions :')
        ion_pool = list(cross_annotations.index)
        court_table = pd.DataFrame(columns = ['ion_list', 'court_size'])
        singletons = list()
        for ion in tqdm(ion_pool):
            if len(cross_annotations.loc[ion].dropna()) == 0:
                singletons.append(ion)
        ion_pool = list(set(ion_pool) - set(singletons))
        total_ions = len(ion_pool)
        while len(ion_pool) > 0 :
            perc = round((1-(len(ion_pool)/total_ions))*100,1)
            sys.stdout.write("\rCourt processing - creating courts : {0}%".format(perc))
            sys.stdout.flush()
            court = [ion_pool[0]]
            court_size = 0
            while court_size != len(court):
                court_size = len(court)
                samples = []
                for ion in court:
                    samples += list(cross_annotations.loc[ion].dropna().index)
                if len(samples) == 0: sys.exit("stuff happenned here")
                samples = list(set(samples))
                new_ions= court.copy()
                for sample in samples:
                    for ion in court:
                        new_ions += Get_cohort_ions_2(ion, sample)
                court = list(set(new_ions))
            court.sort()
            tmp_idx = len(court_table)
            court_table.loc[tmp_idx, 'ion_list'] = court
            court_table.loc[tmp_idx, 'court_size'] = len(court)
            ion_pool = list(set(ion_pool) - set(court))
        court_table.sort_values('court_size', ascending = False, inplace = True)
        court_table.reset_index(drop = True, inplace = True)
        return court_table, singletons
    
    def Cross_neutral_selection(mgf_file):
        """Select the most likeley neutrals from the cross_court_table and report
        the data on the merged node and edge tables.
        """
        modified_cosine = ModifiedCosine(tolerance=mass_error)
        # for each court, find houses:
        print('')
        print("Court processing - selecting houses and reporting results")
        for court in tqdm(cross_court_table.index):
            ion_list = cross_court_table.loc[court, "ion_list"]
            samples = []
            for ion in ion_list:
                samples += list(cross_annotations.loc[ion].dropna().index)
            samples = list(set(samples))
            samples.sort()
            samples_table = pd.DataFrame(index = ion_list, columns = samples)
            for ion in ion_list:
                samples = list(cross_annotations.loc[ion].dropna().index)
                for sample in samples:
                    samples_table.loc[ion, sample] = cross_annotations.loc[ion, sample]
        
            # Produce neutral table
            neutral_table = Cross_Neutral_table(samples_table)
            
            # Do neutral merging:
            merged_neutral_table = Neutral_merging(neutral_table)

            # Eliminate stillborn neutrals:
            drop_list = list()
            for i in merged_neutral_table.index:
                if len(merged_neutral_table.loc[i, adduct_table_merged['Adduct_code']].dropna()) <= 1 :
                    drop_list.append(i)
            merged_neutral_table.drop(drop_list, inplace = True)
            
            ##################################################################
            ################################################## BNR BEGINS HERE
            if run_bnr:
                drop_list = list()
                for i in merged_neutral_table.index:
                    tmp_cos_list = list()
                    adducts = merged_neutral_table.loc[i, adduct_table_merged['Adduct_code']].dropna().index.tolist()
                    if len(bnr_list.intersection(adducts)) == 0:
                        for adduct_1 in adducts:
                            adduct_ions_1 = merged_neutral_table.loc[i, adduct_1]
                            adducts.remove(adduct_1)
                            for ion_1 in adduct_ions_1:
                                ion_1_mgf = merged_node_table.loc[ion_1, "mgf_index"]
                                for adduct_2 in adducts:
                                    adduct_ions_2 = merged_neutral_table.loc[i, adduct_2]
                                    for ion_2 in adduct_ions_2:
                                        ion_2_mgf = merged_node_table.loc[ion_2, "mgf_index"]
                                        score, n_matches = modified_cosine.pair(mgf_file[ion_1_mgf],
                                                                                mgf_file[ion_2_mgf])
                                        if n_matches <= 2 : score = 0.0
                                        tmp_cos_list.append(score)
                        if max(tmp_cos_list) < cosine_hardthreshold:
                            drop_list.append(i)
                merged_neutral_table.drop(drop_list, inplace = True)
                                    
            ################################################## BNR ENDS HERE
            ##################################################################
            
            if len(merged_neutral_table) <= 1 :
                selected_house = list(merged_neutral_table.index)
                merged_neutral_table_selected = merged_neutral_table
            else : 
                end_selection = False
                merged_neutral_table_selected = pd.DataFrame()
                while end_selection == False :
                    # Do ion affinity here:
    ###############################################################################
    ###################################################### ION AFFINITY BEGINS HERE
                    # Ion affinity on neutral remains:
                    # Find conflicted ions
                    contested_ions = pd.DataFrame(columns = ["conflicted_neutrals"])
                    neutral_pool = list(merged_neutral_table.index)
                    while len(neutral_pool) > 1 :
                        neutral_1 = neutral_pool[0]
                        ions_1 = set(flatten(merged_neutral_table.loc[neutral_1, adduct_table_merged['Adduct_code']].dropna()))
                        neutral_pool.remove(neutral_1)
                        for neutral_2 in neutral_pool:
                            ions_2 = set(flatten(merged_neutral_table.loc[neutral_2, adduct_table_merged['Adduct_code']].dropna()))
                            tmp_contested_ions = list(ions_1.intersection(ions_2))
                            for ion in tmp_contested_ions :
                                if ion not in contested_ions.index:
                                    contested_ions.loc[ion, "conflicted_neutrals"] = [neutral_1, neutral_2]
                                    contested_ions.loc[ion, "conflicted_neutrals"].sort()
                                else:
                                    contested_ions.loc[ion, "conflicted_neutrals"] += [neutral_1, neutral_2]
                                    contested_ions.loc[ion, "conflicted_neutrals"] = list(set(contested_ions.loc[ion, "conflicted_neutrals"]))
                                    contested_ions.loc[ion, "conflicted_neutrals"].sort()
                    # Get unresolved neutrals (neutrals with only contested ions)
                    conflicted_neutrals = list()
                    for neutral in merged_neutral_table.index:
                        resolved_ions = set(flatten(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()))
                        resolved_ions = resolved_ions - set(contested_ions.index)
                        if len(resolved_ions) == 0 : conflicted_neutrals.append(neutral)
                        
                    drop_list = list()
                    for ion in contested_ions.index:
                        contested_ions.loc[ion, "conflicted_neutrals"] = list(set(contested_ions.loc[ion, "conflicted_neutrals"]) - set(conflicted_neutrals))
                        if len(contested_ions.loc[ion, "conflicted_neutrals"]) <= 1 :
                            drop_list.append(ion)
                    contested_ions.drop(drop_list, inplace = True)
                    
                    # For each conflicted ion, do ion affinity:
                    same_group_list = list()
                    migration_list = list()
                    resolved_list = list()
                    for ion_1 in contested_ions.index:
                        ion_1_mgf_idx = merged_node_table.loc[ion_1, "mgf_index"]
                        affinity_table = list()
                        resolved = True
                        for neutral in contested_ions.loc[ion_1, "conflicted_neutrals"]:
                            #neutral = contested_ions.loc[ion_1, "conflicted_neutrals"][0]
                            for i in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna().index:
                                if ion_1 in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']][i] : 
                                    ion_1_group = i
                            ion_1_group = adduct_table_merged['Group'][adduct_table_merged['Adduct_code'] == ion_1_group].iloc[0]
                            conflicted_score = False
                            neutral_ions = set(flatten(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()))
                            neutral_ions = list(neutral_ions - set(contested_ions.index))
                            if len(neutral_ions) == 0 :
                                neutral_ions = list(flatten(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()))
                                neutral_ions.remove(ion_1)
                                conflicted_score = True
                            score_table = list()
                            # Get best cosine score with associated n_matches
                            for ion_2 in neutral_ions:
                                for i in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna().index:
                                    if ion_2 in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']][i] : 
                                        ion_2_group = i
                                ion_2_group = adduct_table_merged['Group'][adduct_table_merged['Adduct_code'] == ion_2_group].iloc[0]
                                # cosine, shared peaks, rule_points, dmz, drt, cohort_size
                                ion_2_mgf_idx = merged_node_table.loc[ion_2, "mgf_index"]
                                score, n_matches = modified_cosine.pair(mgf_file[ion_1_mgf_idx],
                                                                        mgf_file[ion_2_mgf_idx])
                                if score < cosine_threshold : score = 0.0
                                prod = score * n_matches
                                if ion_1_group == ion_2_group :
                                    same_group = True
                                else:
                                    same_group = False
                                score_table.append((score, n_matches, prod, same_group))
                            score_table = pd.DataFrame(score_table, columns = ['cos', 'matches', 'prod', 'same_group'])
                            score_table.sort_values("prod", ascending = False, inplace = True)
                            # Get d_mz and d_rt from neutral
                            annotation = merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()
                            ion_count = len(annotation)
                            annotation = [i for i in annotation.index if ion_1 in annotation[i]][0]
                            add_idx = adduct_table_merged.index[adduct_table_merged['Adduct_code'] == annotation][0]
                            ion_neutral = Neutral_mass_calculator(merged_node_table.loc[ion_1, "mz"],
                                                    adduct_table_merged.loc[add_idx, "Adduct_mass"],
                                                    adduct_table_merged.loc[add_idx, "Mol_multiplier"],
                                                    adduct_table_merged.loc[add_idx, "Charge"])
                            d_mz = abs(ion_neutral - merged_neutral_table.loc[neutral, "mass"])
                            d_rt = abs(merged_node_table.loc[ion_1, "rt"] - merged_neutral_table.loc[neutral, "rt"])
                            Species_rule = Validator_choice(annotation, ion_mode)
                            rule_point = Species_rule(ion_1_mgf_idx)
                            if len(score_table) != score_table['same_group'].sum():
                                same_group = False
                            else:
                                same_group = True
                            affinity_table.append((neutral,
                                                   score_table['cos'].iloc[0],
                                                   score_table['matches'].iloc[0],
                                                   rule_point,
                                                   (score_table['cos'].iloc[0] * score_table['matches'].iloc[0]) + rule_point,
                                                   d_mz,
                                                   d_rt,
                                                   d_mz * d_rt,
                                                   ion_count,
                                                   conflicted_score,
                                                   same_group))
                        affinity_table = pd.DataFrame(affinity_table, columns = ['neutral', 'cos',
                                                                                 'matches', 'rule_points', 'score_1',
                                                                                 'd_mz', 'd_rt', 'score_2', 'ion_count',
                                                                                 'conflicted_score', 'same_group'])
                        affinity_table.sort_values(by = ['score_1', 'score_2', 'ion_count'], ascending = [False, True, False], inplace = True)
                        max_score = affinity_table['score_1'].max()
                        if sum(affinity_table['score_1'].between(max_score - 0.05, max_score + 0.05, inclusive = "both")) > 1 :
                            resolved = False
                        migration_list.append(list(affinity_table['neutral']))
                        resolved_list.append(resolved)
                        if len(affinity_table) != affinity_table['same_group'].sum():
                            same_group_list.append(False)
                        else:
                            same_group_list.append(True)
                    contested_ions['selected_neutrals'] = migration_list
                    contested_ions['resolved'] = resolved_list
                    contested_ions['same_group'] = same_group_list
                    contested_ions = contested_ions[contested_ions['same_group']]
                    neutral_ion_count = pd.DataFrame(columns = ['ion_count'])
                    for i in merged_neutral_table.index:
                        neutral_ion_count.loc[i, 'ion_count'] = len(list(flatten(merged_neutral_table.loc[i, adduct_table_merged['Adduct_code']].dropna())))
                    neutral_ion_count["ion_count"] = neutral_ion_count["ion_count"].astype(int)
                    
                    eliminating_neutrals = True
                    while eliminating_neutrals:                
                        neutral_ion_count["new_ion_count"] = neutral_ion_count["ion_count"].copy()
                        for i in contested_ions.index:
                            if contested_ions.loc[i, "resolved"]:
                                loosing_neutrals = contested_ions.loc[i, "selected_neutrals"][1:]
                                for j in loosing_neutrals:
                                    neutral_ion_count.loc[j, "new_ion_count"] = neutral_ion_count.loc[j, "ion_count"] - 1
                        dead_neutrals = neutral_ion_count.index[neutral_ion_count['new_ion_count'] <= 1]
                        if len(dead_neutrals) > 0 :
                            for i in dead_neutrals:
                                neutral_ion_count.drop(i, inplace = True)
                                merged_neutral_table.drop(i, inplace = True)
                                for j in contested_ions.index:
                                    if i in contested_ions.loc[j, "conflicted_neutrals"]:
                                        contested_ions.loc[j, "conflicted_neutrals"].remove(i)
                                    if i in contested_ions.loc[j, "selected_neutrals"] :
                                        contested_ions.loc[j, "selected_neutrals"].remove(i)
                            drop_list = []
                            for i in contested_ions.index:
                                if len(contested_ions.loc[i, "selected_neutrals"]) == 0 :
                                    drop_list.append(i)
                                if len(contested_ions.loc[i, "selected_neutrals"]) == 1 :
                                    contested_ions.loc[i, "resolved"] = True
                            contested_ions.drop(drop_list, inplace = True)
                        else:
                            eliminating_neutrals = False
                    #Once this is done, do the migration:
                    for i in contested_ions.index:
                        for j in contested_ions.loc[i, "selected_neutrals"][1:]:
                            for k in adduct_table_merged['Adduct_code']:
                                if merged_neutral_table.loc[j, k] != None and i in merged_neutral_table.loc[j, k]:
                                    merged_neutral_table.loc[j, k].remove(i)
                                    if len(merged_neutral_table.loc[j, k]) == 0:
                                        merged_neutral_table.loc[j, k] = None
    ######################################################## ION AFFINITY ENDS HERE
    ###############################################################################                
                    
                    # Check neutral compatibility:
                    merged_neutral_table = Neutral_compatibility(merged_neutral_table)
                    
                    # Eject inert neutrals:
                    inert_neutrals = list()
                    neutral_count = len(merged_neutral_table.index)
                    for i in merged_neutral_table.index:
                        inert_neutrals.append((i, len(merged_neutral_table.loc[i, "compatibles"]) + 1))
                    inert_neutrals = pd.DataFrame(inert_neutrals, columns = ['neutral', 'compatible_count'])
                    inert_neutrals = inert_neutrals[inert_neutrals['compatible_count'] == neutral_count]
                    if len(inert_neutrals) > 0 : 
                        merged_neutral_table_selected = merged_neutral_table_selected.append(merged_neutral_table.loc[inert_neutrals['neutral']])
                        merged_neutral_table.drop(inert_neutrals['neutral'], inplace = True)

                    # Check neutral compatibility:
                    merged_neutral_table = Neutral_compatibility(merged_neutral_table)

                    if len(merged_neutral_table) == 0 :
                        selected_house = list(merged_neutral_table_selected.index)
                        end_selection = True
                    else:
                        
                        # Detect uncompatibles (Neutrals compatible with nothing else):
                        uncompatibles = []
                        for i in merged_neutral_table.index:
                            if len(merged_neutral_table.loc[i, "compatibles"]) == 0:
                                uncompatibles.append([i])   
                
                        # Get houses (combinations of compatible neutrals from the same court)
                        houses = []
                        neutral_pool = list(merged_neutral_table.index)
                        tmp_merged_neutral_table = merged_neutral_table.copy()
                        while len(neutral_pool) > 0:
                            neutral_seed = neutral_pool[0]
                            #for neutral_seed in merged_neutral_table.index:
                            next_list, house_list, compatibles_list, left_list = Initialise_house_table_neutrals(neutral_seed, tmp_merged_neutral_table)       
                            if len(house_list) == 0 : houses.append([neutral_seed])
                            while len(house_list) > 0 :
                                next_list, house_list, compatibles_list, left_list, houses = Clean_house_table_neutrals(next_list, house_list, compatibles_list, left_list, houses)
                                next_list, house_list, compatibles_list, left_list = Expand_house_table_neutrals(next_list, house_list, compatibles_list, left_list, tmp_merged_neutral_table)
                            tmp_merged_neutral_table.drop(neutral_seed, inplace = True)
                            for i in tmp_merged_neutral_table.index:
                                if neutral_seed in tmp_merged_neutral_table.loc[i, "compatibles"]:
                                    tmp_merged_neutral_table.loc[i, "compatibles"].remove(neutral_seed)
                            neutral_pool = list(tmp_merged_neutral_table.index)
                        houses += uncompatibles
                        
                        # Create house table for the neutral combinations and assign scores to houses
                        tmp_points_list = list()
                        for i in range(len(houses)):
                            tmp_points = 0.0
                            house = houses[i]
                            for j in house:
                                tmp_points += merged_neutral_table.loc[j, "points"]
                            tmp_points_list.append(tmp_points)
                        house_table = pd.DataFrame()
                        house_table['house'] = houses
                        house_table['points'] = tmp_points_list
                        
                        # Select the best scoring house and report annotation results:
                        selected_house = house_table['points'].idxmax()
                        selected_house = house_table.loc[selected_house, "house"]
                        if len(selected_house) == len(merged_neutral_table.index) :
                            end_selection = True
                        else: 
                            # Ion affinity on neutral remains:
                            # Find conflicted ions
                            contested_ions = pd.DataFrame(columns = ["conflicted_neutrals"])
                            neutral_pool = list(merged_neutral_table.index)
                            while len(neutral_pool) > 1 :
                                neutral_1 = neutral_pool[0]
                                ions_1 = set(flatten(merged_neutral_table.loc[neutral_1, adduct_table_merged['Adduct_code']].dropna()))
                                neutral_pool.remove(neutral_1)
                                for neutral_2 in neutral_pool:
                                    ions_2 = set(flatten(merged_neutral_table.loc[neutral_2, adduct_table_merged['Adduct_code']].dropna()))
                                    tmp_contested_ions = list(ions_1.intersection(ions_2))
                                    for ion in tmp_contested_ions :
                                        if ion not in contested_ions.index:
                                            contested_ions.loc[ion, "conflicted_neutrals"] = [neutral_1, neutral_2]
                                            contested_ions.loc[ion, "conflicted_neutrals"].sort()
                                        else:
                                            contested_ions.loc[ion, "conflicted_neutrals"] += [neutral_1, neutral_2]
                                            contested_ions.loc[ion, "conflicted_neutrals"] = list(set(contested_ions.loc[ion, "conflicted_neutrals"]))
                                            contested_ions.loc[ion, "conflicted_neutrals"].sort()
                            
                            winning_neutral = list()
                            loosing_neutrals = list()
                            for i in contested_ions.index:
                                tmp_win = list(set(selected_house).intersection(contested_ions.loc[i, "conflicted_neutrals"]))
                                tmp_loss = list(set(contested_ions.loc[i, "conflicted_neutrals"]) - set(tmp_win))
                                winning_neutral.append(tmp_win)
                                loosing_neutrals.append(tmp_loss)
                            contested_ions['winning_neutral'] = winning_neutral
                            contested_ions['loosing_neutrals'] = loosing_neutrals
                            for i in contested_ions.index:
                                for neutral in contested_ions.loc[i, "loosing_neutrals"]:
                                    for adduct in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna().index:
                                        merged_neutral_table.loc[neutral, adduct] = list(set(merged_neutral_table.loc[neutral, adduct]) - set([i]))
                                        if len(merged_neutral_table.loc[neutral, adduct]) == 0 :
                                            merged_neutral_table.loc[neutral, adduct] = None
                            dead_neutrals = list()
                            for neutral in merged_neutral_table.index:
                                if len(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()) <= 1 :
                                    dead_neutrals.append(neutral)
                            merged_neutral_table.drop(dead_neutrals, inplace = True)

                                       
            # Report results
            merged_neutral_table = merged_neutral_table_selected
            for neutral in selected_house:
                neutral_mass = merged_neutral_table.loc[neutral, "mass"]
                neutral_rt = merged_neutral_table.loc[neutral, "rt"]
                neutral_idx = merged_node_table.index.max() + 1
                merged_node_table.loc[neutral_idx] = [neutral_mass, neutral_rt, 0, 0, None, "neutral", None]
                annotations = merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()
                for annotation in annotations.index:
                    ion_list = annotations[annotation]
                    for ion in ion_list:
                        mz_gap = round(merged_node_table.loc[ion, "mz"] - neutral_mass, 4)
                        new_edge = merged_edge_table.index.max() + 1
                        merged_edge_table.loc[new_edge] = [neutral_idx, ion, 0, 0, 0.0, 0.0, mz_gap, "add_edge", None, annotation, annotation]
                        merged_node_table.loc[neutral_idx, 'TIC'] += merged_node_table.loc[ion, 'TIC']
                        if merged_node_table.loc[ion, 'Adnotation'] == None :
                            merged_node_table.loc[ion, 'Adnotation'] = annotation
                        else :
                            sys.exit("ERROR HERE FOR SOME REASON : already adnotated ion?")
        return
    
    def Cross_Neutral_table(samples_table):
        """Produces a neutral_table, containing among other things, the computed
        neutral mass for each ion in the samples table.
        """
        neutral_table = pd.DataFrame(columns = ['ion_ID', 'mz', 'rt', 'annotation', 'neutral'])
        for ion in samples_table.index:
            annotations = samples_table.loc[ion].dropna().unique()
            mz = merged_node_table.loc[ion, "mz"]
            rt = merged_node_table.loc[ion, "rt"]
            for annotation in annotations :
                tmp_sample = cross_annotations.columns[cross_annotations.loc[ion] == annotation][0]
                neutral = float(cross_neutrals.loc[ion, tmp_sample])
                neutral_table.loc[len(neutral_table)] = [ion, mz, rt, annotation, neutral]
        return neutral_table
    
    def Neutral_merging(neutral_table):
        """Merges neutrals from the neutral table based on their retention time
        and molecular masses. Returns a merged neutral table, with the neutral's
        retention time, points accumulated, the ions involved and their annotations.
        """
        neutral_table_2 = pd.DataFrame(columns = ["mass", "rt", "points"] + list(adduct_table_merged['Adduct_code']))
        neutral_idx = neutral_table.index
        while len(neutral_idx) > 1:
            neutral_pool = [neutral_idx[0]]
            pool_size = 0
            pool_table = neutral_table.loc[neutral_pool]
            while len(neutral_pool) != pool_size:
                pool_size = len(neutral_pool)
                tmp_rt = (pool_table['rt'].min() + pool_table['rt'].max())/2
                tmp_neutral = (pool_table['neutral'].min() + pool_table['neutral'].max())/2
                pool_table = neutral_table[neutral_table['rt'].between(tmp_rt - rt_error, tmp_rt + rt_error)]
                pool_table = pool_table[pool_table['neutral'].between(tmp_neutral - mass_error, tmp_neutral + mass_error)]
                neutral_pool = list(pool_table.index)
            tmp_idx = len(neutral_table_2)
            neutral_table_2.loc[tmp_idx] = [None]*len(neutral_table_2.columns)
            neutral_table_2.loc[tmp_idx, "mass"] = round(pool_table['neutral'].mean(), 4)
            neutral_table_2.loc[tmp_idx, "rt"] = round(pool_table['rt'].mean(), 2)
            for annotation in pool_table['annotation'].unique():
                neutral_table_2.loc[tmp_idx, annotation] = list(pool_table['ion_ID'][pool_table["annotation"] == annotation])
            adducts = list(neutral_table_2.loc[tmp_idx, adduct_table_merged['Adduct_code']].dropna().index)
            tmp_points = 0
            for adduct in adducts:
                tmp_points += 1/(adduct_table_merged["Complexity"][adduct_table_merged["Adduct_code"] == adduct].iloc[0])
            neutral_table_2.loc[tmp_idx, "points"] = tmp_points
            neutral_idx = list(set(neutral_idx) - set(neutral_pool))
        return neutral_table_2   
    
    def Neutral_compatibility(merged_neutral_table):
        """For each neutral in the neutral table, adds a list of compatible neutrals
        found in the same table.
        """
        merged_neutral_table['compatibles'] = [None]*len(merged_neutral_table)
        for i in merged_neutral_table.index:
            merged_neutral_table.at[i, "compatibles"] = list()
        for neutral_1 in merged_neutral_table.index:
            annotations_1 = merged_neutral_table.loc[neutral_1, adduct_table_merged['Adduct_code']].dropna()
            annotations_1 = set(flatten(annotations_1))      
            for neutral_2 in merged_neutral_table.index:
                if neutral_1 == neutral_2 : continue
                annotations_2 = merged_neutral_table.loc[neutral_2, adduct_table_merged['Adduct_code']].dropna()
                annotations_2 = set(flatten(annotations_2))
                if len(annotations_1.intersection(annotations_2)) == 0 :
                    merged_neutral_table.loc[neutral_2, 'compatibles'].append(neutral_1)
        for i in merged_neutral_table.index:
            merged_neutral_table.loc[i, "compatibles"] = list(set(merged_neutral_table.loc[i, "compatibles"]))
            merged_neutral_table.loc[i, "compatibles"].sort()        
        return merged_neutral_table
    
    def Initialise_house_table_neutrals(neutral_seed, merged_neutral_table):
        """Starts house building by combining compatible cohorts, using an initial
        cohort as a neutral_seed.
        """
        next_list = merged_neutral_table.loc[neutral_seed, "compatibles"].copy()
        house_list = [[neutral_seed]]*len(next_list)
        compatibles_list = [None]*len(next_list)
        left_list = [0]*len(next_list)
        added = []
        for i in range(len(next_list)):
            added.append(next_list[i])
            tmp_compatibles = list(set(next_list).intersection(merged_neutral_table.loc[next_list[i], "compatibles"]))
            tmp_compatibles = list(set(tmp_compatibles) - set(added))
            compatibles_list[i] = tmp_compatibles
            left_list[i] = len(compatibles_list[i])
        return next_list, house_list, compatibles_list, left_list
        
        # Clean house table and export results
    def Clean_house_table_neutrals(next_list, house_list, compatibles_list, left_list, houses):
        """Removes houses that finished expanding (no compatible cohorts left) from
        the house table, and exports them to the houses list.
        """
        export_idx = [i for i in range(len(left_list)) if left_list[i] == 0]
        export_idx.sort(reverse = True)
        for i in export_idx:
            new_house = house_list[i] + [next_list[i]]
            new_house = list(set(new_house))
            new_house.sort()
            houses.append(new_house) 
            next_list.pop(i)
            house_list.pop(i)
            compatibles_list.pop(i)
            left_list.pop(i)
        return next_list, house_list, compatibles_list, left_list, houses
    
    def Expand_house_table_neutrals(next_list, house_list, compatibles_list, left_list, merged_neutral_table):
        """Expands the houses in the house table by appending to them one of each
        possible compatible cohorts.
        """
        drop_list = []
        for i in range(len(next_list)):
            drop_list.append(i)
            new_house = house_list[i] + [next_list[i]]
            added = []
            compatibles = compatibles_list[i]
            for j in compatibles:
                #j = compatibles[0]
                added.append(j)
                next_list.append(j)
                house_list.append(new_house)
                tmp_compatibles = list(set(compatibles).intersection(merged_neutral_table.loc[j, "compatibles"]))
                tmp_compatibles = list(set(tmp_compatibles) - set(added))
                compatibles_list.append(tmp_compatibles)
                left_list.append(len(set(compatibles).intersection(merged_neutral_table.loc[j, "compatibles"])))
        drop_list.sort(reverse = True)
        for i in drop_list:
            next_list.pop(i)
            house_list.pop(i)
            compatibles_list.pop(i)
            left_list.pop(i)
        return next_list, house_list, compatibles_list, left_list  
    
    def Update_node_table(adduct_table_merged):
        """Updates the node table with the "adduct" status on adnotated nodes.
        Also adds an "adduct_count" columns with the number of ions directly
        connected to the neutral node.
        """
        print("Adding rule points...")
        rule_points = list()
        for i in tqdm(merged_node_table.index):
            adnotation = merged_node_table.loc[i, 'Adnotation']
            if adnotation == None : 
                rule_points.append(0)
            else:
                sample = (cross_annotations.loc[i].dropna() == adnotation).index[0]
                rule_points.append(int(cross_rules.loc[i, sample]))
        merged_node_table.insert(merged_node_table.columns.get_loc("mgf_index") + 1, 'rule_points', rule_points)

        
        adnoted_ions = merged_node_table['Adnotation'].dropna().index
        print("Updating node table status")
        for i in tqdm(adnoted_ions):
            merged_node_table.loc[i, "status"] = "adduct"
        
        merged_node_table['adduct_count'] = [0]*len(merged_node_table)
        for i in merged_node_table.index[merged_node_table['status'] == "neutral"]:
            count = sum(merged_edge_table['node_1'] == i)
            merged_node_table.loc[i, "adduct_count"] = count
        
        for adduct in tqdm(merged_node_table['Adnotation'].dropna().unique()):
            clean_adduct = adduct_table_merged['Adduct'][adduct_table_merged['Adduct_code'] == adduct].iloc[0]
            change_idx = merged_node_table.index[merged_node_table['Adnotation'] == adduct]
            for i in change_idx : 
                merged_node_table.loc[i, "Adnotation"] = clean_adduct
        merged_node_table['mz'] = merged_node_table['mz'].round(4)
        merged_node_table['rt'] = merged_node_table['rt'].round(3)
        adnotation = merged_node_table['Adnotation'].copy()
        merged_node_table.drop('Adnotation', axis = 1, inplace = True)
        merged_node_table['Adnotation'] = adnotation
        return
    
    def Update_edge_table(merged_edge_table, adduct_table_merged):
        """Updates edge table with singleton nodes
        """
        print("Updating edge table...")
        singleton_ions = merged_node_table.index[merged_node_table["status"] == "singleton"]
        singleton_df = pd.DataFrame()
        singleton_df['node_1'] = list(singleton_ions)
        singleton_df['node_2'] = list(singleton_ions)
        singleton_df['matched_peaks'] = [0]*len(singleton_ions)
        singleton_df['total_peaks'] = [0]*len(singleton_ions)
        singleton_df['matching_score'] = [0.0]*len(singleton_ions)
        singleton_df['rt_gap'] = [0.0]*len(singleton_ions)
        singleton_df['mz_gap'] = [0.0]*len(singleton_ions)
        singleton_df['status'] = ["self_edge"]*len(singleton_ions)
        singleton_df['Fragnotation'] = [None]*len(singleton_ions)
        singleton_df['All_annotations'] = [None]*len(singleton_ions)
        singleton_df['Adnotation'] = [None]*len(singleton_ions)
        merged_edge_table = merged_edge_table.append(singleton_df, ignore_index = True)
        print("Replacing adduct codes by adduct formulas...")
        for adduct in tqdm(merged_edge_table['Adnotation'].dropna().unique()):
            clean_adduct = adduct_table_merged['Adduct'][adduct_table_merged['Adduct_code'] == adduct].iloc[0]
            change_idx = merged_edge_table.index[merged_edge_table['Adnotation'] == adduct]
            for i in change_idx : 
                merged_edge_table.loc[i, "Adnotation"] = clean_adduct
                merged_edge_table.loc[i, "All_annotations"] = clean_adduct
        merged_edge_table['mz_gap'] = merged_edge_table['mz_gap'].astype(float)
        merged_edge_table['mz_gap'] = merged_edge_table['mz_gap'].round(4)
        merged_edge_table['rt_gap'] = merged_edge_table['rt_gap'].round(3)
        all_annotations = merged_edge_table["All_annotations"].copy()
        merged_edge_table.drop('All_annotations', axis = 1, inplace = True)
        merged_edge_table['All_annotations'] = all_annotations
        return merged_edge_table
    
    def Samplewise_export(full_csv_file, out_path_samples, merged_edge_table, merged_node_table) : 
        print("Exporting sample-wise tables...")
        full_csv = pd.read_csv(full_csv_file, index_col ="row ID")
        samples = full_csv.columns
        samples = samples.str.replace(mzmine_suffix, "", regex = False)
        full_csv.columns = samples
        samples = list(samples.drop(["row m/z", "row retention time"]))
        
        if not os.path.isdir(out_path_samples):
            os.mkdir(out_path_samples)
        
        for sample in tqdm(samples):
            # get sample ions
            ions = full_csv.index[full_csv[sample] > 0.0]
        
            # Get sample neutrals
            kept_edges = list()
            neutral_edges = merged_edge_table[merged_edge_table['status'] == "add_edge"]
            for i in neutral_edges.index:
                if neutral_edges.loc[i, "node_2"] in ions :
                    kept_edges.append(i)
        
            # Get ion edges
            ion_edges = merged_edge_table[merged_edge_table['status'] != "add_edge"]
            for i in ion_edges.index:
                if ion_edges.loc[i, "node_1"] in ions:
                    if ion_edges.loc[i, "node_2"] in ions:
                        kept_edges.append(i)
            kept_edges.sort()
            sample_edges = merged_edge_table.loc[kept_edges]
            sample_edges.sort_values('node_1', inplace = True)
            sample_edges.reset_index(inplace = True, drop = True)
            
            kept_nodes = list(set(list(sample_edges['node_1']) + list(sample_edges['node_2'])))
            kept_nodes.sort()
            sample_nodes = merged_node_table.loc[kept_nodes]
            
            sample_nodes.to_csv(out_path_samples + sample+"_nodes.csv", index_label = "feature_index")
            sample_edges.to_csv(out_path_samples + sample+"_edges.csv", index_label = "Index")
        return
    
    ###############################################################################
    ############################## START ADNOTATION ###############################
    ###############################################################################

    import os
    import pandas as pd
    import numpy as np
    import sys
    from tqdm import tqdm
    from pandas.core.common import flatten
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.similarity import ModifiedCosine
    
    # Load parameters:
    mass_error= params['an_mass_error']
    prec_mass_error= params['an_prec_mass_error']
    rt_error= params['an_rt_error']
    cosine_threshold= params['an_cos_threshold']      
    cosine_hardthreshold=  params['an_hardcos_threshold']
    run_bnr = params['an_run_bnr']
    mzmine_suffix= params['mzmine_suffix']
    if ion_mode == "NEG":
        in_path_full= params['neg_out_0']
        csv_file= params['neg_csv']
        mgf_file= params['neg_mgf']
        in_path_mgf= params['neg_out_1']
        in_path_csv= params['neg_out_2']
        out_path_full= params['neg_out_3_1']
        out_path_samples= params['neg_out_3_2']
        adduct_table_primary= params['an_addtable_primary_neg']
        adduct_table_secondary= params['an_addtable_secondary_neg']
        bnr_list = set(params['an_bnr_neg'])
    elif ion_mode == "POS":
        in_path_full= params['pos_out_0']
        csv_file= params['pos_csv']
        mgf_file= params['pos_mgf']
        in_path_mgf= params['pos_out_1']
        in_path_csv= params['pos_out_2']
        out_path_full= params['pos_out_3_1']
        out_path_samples= params['pos_out_3_2']
        adduct_table_primary= params['an_addtable_primary_pos']
        adduct_table_secondary= params['an_addtable_secondary_pos']
        bnr_list = set(params['an_bnr_pos'])
    else:
        print('Ion mode must be either "NEG" or "POS"')
        return
    
    full_csv_file = in_path_full + csv_file
    full_mgf_file = in_path_full + mgf_file
    
    input_files = pd.DataFrame()
    input_files['mgf_file'] = os.listdir(in_path_mgf)
    input_files['base_name'] = input_files['mgf_file'].copy().str.replace(f'{ion_mode}_', "", regex = False)
    input_files['base_name'] = input_files['base_name'].str.replace('.mgf', '', regex = False)
    
    adduct_table_primary = pd.read_csv("./params/" + adduct_table_primary, sep = "\t")
    adduct_table_secondary = pd.read_csv("./params/" + adduct_table_secondary, sep = "\t")
    adduct_table_merged= adduct_table_primary.append(adduct_table_secondary, ignore_index = True)
    
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    
    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
    
    if os.path.isfile(out_path_full + "cross_sample_annotations.csv"):
        cross_annotations = pd.read_csv(out_path_full + "cross_sample_annotations.csv", index_col = "Feature_ID", dtype = str)
        cross_annotations = cross_annotations.replace({np.nan: None})
        cross_courts = pd.read_csv(out_path_full + "cross_sample_courts.csv", index_col = "Feature_ID")
        cross_courts = cross_courts.replace({np.nan: None})
        cross_houses = pd.read_csv(out_path_full + "cross_sample_houses.csv", index_col = "Feature_ID", dtype = str)
        cross_houses = cross_houses.replace({np.nan: None})
        cross_points = pd.read_csv(out_path_full + "cross_sample_points.csv", index_col = "Feature_ID", dtype = float)
        cross_points = cross_points.replace({np.nan: None})
        cross_rules = pd.read_csv(out_path_full + "cross_sample_rules.csv", index_col = "Feature_ID", dtype = str)
        cross_rules = cross_rules.replace({np.nan: None})
        cross_neutrals = pd.read_csv(out_path_full + "cross_sample_neutrals.csv", index_col = "Feature_ID", dtype = str)
        cross_neutrals = cross_neutrals.replace({np.nan: None})
        feature_ids = Get_feature_ids(full_mgf_file)
    else:
        cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals, feature_ids = Cross_sample_tables(full_mgf_file)
    
    for x in input_files.index :
        sample_base_name = input_files.loc[x, "base_name"]
        if sample_base_name in cross_annotations.columns : continue
        
        file_basename = input_files.at[x, 'base_name']
        print("Processing " + file_basename)
        
        edge_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_edges.csv", 
                                 index_col = 'Index')
        node_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_nodes.csv", 
                                index_col = 'Index')
        full_merged_table = pd.DataFrame()
        
        mgf_file = in_path_mgf + input_files.loc[x, "mgf_file"]
        mgf_file = list(load_from_mgf(mgf_file))
        mgf_file = [Spectrum_processing(s) for s in mgf_file]
        
        duplicate_df = pd.DataFrame() # Stores duplicate spectra
    
        for i in tqdm(node_table.index) :
            # Load data for first ion to be analysed
            #if node_table.loc[i, 'status'] == "fragment" : continue
            if abs(int(node_table.loc[i, 'charge'])) > 1 : continue
            ion1_rt = node_table.loc[i, 'rt']
            ion1_mz = node_table.loc[i, 'mz']
            #ion1_id = node_table.loc[i, 'feature_id']
    
            # Slice the DF to keep only ions coeluting with current ion
            coelution_table = Rt_slicer(ion1_rt, rt_error, i, node_table)
            
            # Filter the table to remove multicharged ions
            coelution_table = Charge_filter(1, coelution_table)
        
            # Considering all adduct possibilities for current ion, 
            #calculate the associated hypothetical molecular mass:
            neutral_table = Neutral_table(ion1_mz, adduct_table_primary)
        
            # Considering the first ionisation hypotheses and the associated 
            # neutrals, all other possible ions are calculated based on the 
            # hypothetical molecule. These are stored in an ion_hypotheses_table:
            # Ion1_index: the first ionisation hypothesis from the previous 
            # neutral_table's index, 
            # Ion2_index: the second ionisation hypothesis considering the first
            # hypothesis and the hypothetical molecule.
            # Ion2_complexity: Complexity score from the Ion2, used for later calculations
            # Ion2_mz: Ion2 m/z value, which is to be searched in the coelution
            # table to validate the first ionisation hypothesis.
            ion_hypotheses_table = Ion_hypotheses_table(neutral_table)
        
        
            # Award points for each combination of hypotheses (ionisation 1, neutral
            # 1 and ionisation 2) based on the existance of the Ion2 in the 
            # coelution table.
            ion_hypotheses_table = Point_counter(ion_hypotheses_table, mass_error)
    
        
            # the ion_hypotheses_table is filtered to keep only level 2 hypotheses that were detected
            # Also, when multiple hits are found in one single RT slice, their numbers are 
            # reported in the 'hits' column. However it wouldn't make sense for all these
            # hits to have the same m/z value and be different level 2 adducts for the same
            # level 1 molecule, so there should be at most one of these m/z values that is correct
            # hence, the maximum number of hits is set to 1
            #ion_hypotheses_table = Hit_filter(ion_hypotheses_table)
        
            # Cosine validation: the cosine score will count as additional points
            # unless it is below cosine threshold, in which case the score will be 
            # reduced to 0 and the hypothesis will be eliminated if the ion pairs
            # belong to the same adduct group (H, Na, K, Cl)
            ion_hypotheses_table, tmp_duplicates = Cosine_validation(i, cosine_threshold, mgf_file)
            duplicate_df = duplicate_df.append(tmp_duplicates, ignore_index = True)
        
            # Transfer the number of hits to each level 1 ion hypothesis in level_1_ion_hypothesis_table 
            # Points are weighted in the following fashion : (number of hits / level1_complexity)*defragHit*l2frag
            #What is defraghit? what is l2frag?
            # One thing might be checking if the current level 1 ion hypothesis is a fragment from another ion 
            # that should be present in the coeluted ions
            
            # Species rules points : award points for annotations that can be confirmed
            ion_hypotheses_table = Species_rules(i, ion_hypotheses_table)
    
            
        
            
            # Solvent complex confirmation :
            # award points to level 2 hypotheses with solvent complexes for which the
            # uncomplexed form is present in the coeluted table. The uncomplexed form
            # is also rewarded
            ion_hypotheses_table = Complex_points(neutral_table, ion_hypotheses_table)
            
    
            # In case of multiple hits for ion_2, the ion with the closest m/z value
            # to the target is chosen.
            #ion_hypotheses_table = Multi_hit_correcter(ion_hypotheses_table, coelution_table)
            
            # Do a fragment confirmation if the pairs are already linked as fragments
            # (more chances of them being actually related adducts)
            
            ion_hypotheses_table = Fragnotator_points(i, ion_hypotheses_table,
                                                              edge_table, node_table.loc[i, 'status'])
    
            # Get the adduct code for each ionisation hypothesis
            ion_hypotheses_table = Get_adduct_code(ion_hypotheses_table, neutral_table)
            
    
            # Now, count points using hits, solvent points, fragment points and complexity
            ion_hypotheses_table = Weighted_points(ion_hypotheses_table)
        
            
            # Convert the level 2 table into a true adduct identity probability table,
            merged_adducts_table = Merged_adducts_table(ion_hypotheses_table, i)
    
            if len(merged_adducts_table.index) == 0 : continue
            full_merged_table = full_merged_table.append(merged_adducts_table)
            full_merged_table.reset_index(inplace = True, drop = True)
        
        cohort_table = Cohort_Table(neutral_table, full_merged_table, node_table)
        
        supercohorts_table, transition_table = Supercohorts_table(cohort_table)
        
        court_table = Get_Court_Table(supercohorts_table)
        
        court_table = House_selection(court_table)
    
        Cross_sample_report(court_table)
        
        Secondary_adducts(mgf_file)
        
        cross_annotations.to_csv(out_path_full+"cross_sample_annotations.csv", index_label = "Feature_ID")
        cross_courts.to_csv(out_path_full+"cross_sample_courts.csv", index_label = "Feature_ID")
        cross_houses.to_csv(out_path_full+"cross_sample_houses.csv", index_label = "Feature_ID")
        cross_points.to_csv(out_path_full+"cross_sample_points.csv", index_label = "Feature_ID")
        cross_rules.to_csv(out_path_full+"cross_sample_rules.csv", index_label = "Feature_ID")
        cross_neutrals.to_csv(out_path_full+"cross_sample_neutrals.csv", index_label = "Feature_ID")
    
    # Load and merge fragnoted 
    merged_node_table, merged_edge_table = Merged_tables(input_files, feature_ids)
    merged_edge_table['Adnotation'] = [None]*len(merged_edge_table)
    merged_node_table['Adnotation'] = [None]*len(merged_node_table)
    
    # First, resolve the yet sample-specific unresolved ions by cross examining between
    # all samples:
    Annotation_resolver()
    
    # Now rebuild houses at a cross-sample level to check compatibilities between
    # annotations
    cross_court_table, singletons = Cross_court_table()    
    
    # Select the most likeley neutrals from the cross_court_table and report
    # the data on the merged node and edge tables.
    mgf_file = list(load_from_mgf(full_mgf_file))
    mgf_file = [Spectrum_processing(s) for s in mgf_file]
    Cross_neutral_selection(mgf_file)
    
    # Update node table status:
    Update_node_table(adduct_table_merged)
    
    # Update edge table with singletons:
    merged_edge_table = Update_edge_table(merged_edge_table, adduct_table_merged)
    
    # Export merged tables here
    merged_edge_table.to_csv(out_path_full + ion_mode + "_full_edges.csv", index_label = "Index")
    merged_node_table.to_csv(out_path_full + ion_mode + "_full_nodes.csv", index_label = "feature_id")
    
    # Also export sample wise slices of the tables for viewing at this stage
    if params['an_export_samples'] : 
        Samplewise_export(full_csv_file, out_path_samples, merged_edge_table, merged_node_table)
    return
