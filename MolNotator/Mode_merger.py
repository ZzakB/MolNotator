def Moder_merger(params : dict):

    def Solo_M1mHpC4H11N(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 72.081324
        mz_Cl = 34.968853 + mz - 72.081324
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl])
    
    def Solo_M1mHpHCOOH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 44.997654
        mz_Cl = 34.968853 + mz - 44.997654
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl])
    
    def Solo_M1m2HpNapHCOOH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 66.979600
        mz_Cl = 34.968853 + mz - 66.979600
        mz_m2HpNa = 20.97412 + mz - 66.979600
        
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])


    def Solo_M1m2HpNa(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 66.979600
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H])
    
    def Solo_M1m2HpK(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + mz - 36.948058
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H])
    
    def Solo_M2mHpC4H11N(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 72.081324)/2
        mz_Cl = 34.968853 + (mz - 72.081324)/2
        mz_m2HpNa = 20.97412 + (mz - 72.081324)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])
    
    def Solo_M2mHpHCOOH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 44.997654)/2
        mz_Cl = 34.968853 + (mz - 44.997654)/2
        mz_m2HpNa = 20.97412 + (mz - 44.997654)/2
        mz_mHpHCOOH = 44.997654 + (mz - 44.997654)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_mHpHCOOH = peaks.between(mz_mHpHCOOH - prec_mass_error, mz_mHpHCOOH + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_mHpHCOOH])
    
    def Solo_M2mH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz + 1.007825)/2
        mz_Cl = 34.968853 + (mz + 1.007825)/2
        mz_m2HpNa = 20.97412 + (mz + 1.007825)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])
    
    def Solo_M2pCl(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 34.968853)/2
        mz_Cl = 34.968853 + (mz - 34.968853)/2
        mz_m2HpNa = 20.97412 + (mz - 34.968853)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa])
    
    def Solo_M2m2HpNapHCOOH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 66.979600)/2
        mz_Cl = 34.968853 + (mz - 66.979600)/2
        mz_m2HpNa = 20.97412 + (mz - 66.979600)/2
        mz_m2HpNapHCOOH = 66.9796 + (mz - 66.979600)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH])
    
    def Solo_M2m2HpNa(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 20.97412)/2
        mz_Cl = 34.968853 + (mz - 20.97412)/2
        mz_m2HpNa = 20.97412 + (mz - 20.97412)/2
        mz_m2HpNapHCOOH = 66.9796 + (mz - 20.97412)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH])
    
    def Solo_M2m2HpK(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = -1.007825 + (mz - 36.948058)/2
        mz_Cl = 34.968853 + (mz - 36.948058)/2
        mz_m2HpNa = 20.97412 + (mz - 36.948058)/2
        mz_m2HpNapHCOOH = 66.9796 + (mz - 36.948058)/2
        mz_m2HpK = 36.948058 + (mz - 36.948058)/2
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK])
    
    def Solo_M3mH(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])
    
    def Solo_M3pCl(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])
    
    def Solo_M3m2HpNapHCOOH(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])
    
    def Solo_M3m2HpNa(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK])
    
    def Solo_M4mH(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                    valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])
    
    def Solo_M4pCl(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                    valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])
    
    def Solo_M4m2HpNapHCOOH(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                    valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])
    
    def Solo_M4m2HpNa(ion_idx, mgf_file):
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
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Cl = peaks.between(mz_Cl - prec_mass_error, mz_Cl + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNa = peaks.between(mz_m2HpNa - prec_mass_error, mz_m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpNapHCOOH = peaks.between(mz_m2HpNapHCOOH - prec_mass_error, mz_m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_m2HpK = peaks.between(mz_m2HpK - prec_mass_error, mz_m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2mH = peaks.between(mz_M2mH - prec_mass_error, mz_M2mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pCl = peaks.between(mz_M2pCl - prec_mass_error, mz_M2pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNa = peaks.between(mz_M2m2HpNa - prec_mass_error, mz_M2m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpNapHCOOH = peaks.between(mz_M2m2HpNapHCOOH - prec_mass_error, mz_M2m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2m2HpK = peaks.between(mz_M2m2HpK - prec_mass_error, mz_M2m2HpK + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3mH = peaks.between(mz_M3mH - prec_mass_error, mz_M3mH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pCl = peaks.between(mz_M3pCl - prec_mass_error, mz_M3pCl + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNa = peaks.between(mz_M3m2HpNa - prec_mass_error, mz_M3m2HpNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpNapHCOOH = peaks.between(mz_M3m2HpNapHCOOH - prec_mass_error, mz_M3m2HpNapHCOOH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3m2HpK = peaks.between(mz_M3m2HpK - prec_mass_error, mz_M3m2HpK + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Cl, valid_m2HpNa, valid_m2HpNapHCOOH, valid_m2HpK,
                    valid_M2mH, valid_M2pCl, valid_M2m2HpNa, valid_M2m2HpNapHCOOH, valid_M2m2HpK,
                    valid_M3mH, valid_M3pCl, valid_M3m2HpNa, valid_M3m2HpNapHCOOH, valid_M3m2HpK])
    
    def Solo_M2pH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 1.007825)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H])
    
    def Solo_M2pHpCH3CN(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        
        return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])
    
    def Solo_M2pHpCH3OH(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        
        return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])
    
    
    def Solo_M2pHpHCOOH(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])
    
    def Solo_M2pNH4(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH])
    
        
    def Solo_M2pNa(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 22.98977)/2
        mz_Na = 22.98977 + (mz - 22.98977)/2
        mz_NapCH3CN = 64.016319 + (mz - 22.98977)/2
        mz_NapCH3OH = 55.015985 + (mz - 22.98977)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])
    
    def Solo_M2pNapCH3OH(ion_idx, mgf_file) : 
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 55.015985)/2
        mz_Na = 22.98977 + (mz - 55.015985)/2
        mz_NapCH3CN = 64.016319 + (mz - 55.015985)/2
        mz_NapCH3OH = 55.015985 + (mz - 55.015985)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])
    
    def Solo_M2pNapCH3CN(ion_idx, mgf_file) : 
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 64.016319)/2
        mz_Na = 22.98977 + (mz - 64.016319)/2
        mz_NapCH3CN = 64.016319 + (mz - 64.016319)/2
        mz_NapCH3OH = 55.015985 + (mz - 64.016319)/2
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3CN, valid_NapCH3OH])
    
    def Solo_M2pK(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 38.963708)/2
        mz_Na = 22.98977 + (mz - 38.963708)/2
        mz_K = 38.963708 + (mz - 38.963708)/2
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Na, valid_K])
    
    def Solo_M1pHpCH3CN(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 42.034374
        mz_Na = 22.98977 + mz - 42.034374
        mz_HpCH3OH = 33.034040 + mz - 42.034374
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Na, valid_HpCH3OH])
    
    def Solo_M1pHpCH3OH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 33.034040
        mz_Na = 22.98977 + mz - 33.034040
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Na])
    
    def Solo_M1pHpHCOOH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 47.013304
        mz_Na = 22.98977 + mz - 47.013304
        mz_HpCH3OH = 33.034040 + mz - 47.013304
        mz_HpCH3CN = 42.034374 + mz - 47.013304
     
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
    
        return sum([valid_H, valid_Na, valid_HpCH3OH, valid_HpCH3CN])
    
    def Solo_M1pNa(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 22.989770
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum()
    
        return valid_H
    
    def Solo_M1pNapCH3CN(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = mz - 64.016319 + 1.007825
        mz_Na = mz - 64.016319 + 22.98977
        mz_NapCH3OH = mz - 64.016319 + 55.015985
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Na, valid_NapCH3OH])
    
    def Solo_M1pNapCH3OH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 55.015985
        mz_Na = 22.98977 + mz - 55.015985
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_Na])
    
    def Solo_M1pNH4(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 18.034374
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum()
        return valid_H
    
    def Solo_M1pNH4pCH3CN(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 59.060923
        mz_NH4 = 18.034374 + mz - 59.060923
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        return sum([valid_H, valid_NH4])
    
    def Solo_M1pK(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 38.963708
        mz_NH4 =  18.034374 + mz - 38.963708
        mz_Na = 22.98977 + mz - 38.963708
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum()
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum()
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum()
        return sum([valid_H, valid_NH4, valid_Na])
    
    def Solo_M1pKpCH3OH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + mz - 70.989923
        mz_NH4 =  18.034374 + mz - 70.989923
        mz_Na = 22.98977 + mz - 70.989923
        mz_K = 38.963708 + mz - 70.989923
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum()
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum()
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum()
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum()
        return sum([valid_H, valid_NH4, valid_Na, valid_K])
    
    def Solo_M3pH(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 1.007825)/3
        mz_M2pH = 1.007825 + (mz - 1.007825)*(2/3)
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
 
        return sum([valid_H, valid_M2pH])
    
    def Solo_M3pNH4(ion_idx, mgf_file):
        mz = mgf_file[ion_idx].get('pepmass')[0]
        peaks = pd.Series(mgf_file[ion_idx].peaks.mz)
        mz_H = 1.007825 + (mz - 18.034374)/3
        mz_NH4 = 18.034374 + (mz - 18.034374)/3
        mz_M2pH = 1.007825 + (mz - 18.034374)*(2/3)
        mz_M2pNH4 = 18.034374 + (mz - 18.034374)*(2/3)

    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0  
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0         
        return sum([valid_H, valid_NH4, valid_M2pH, valid_M2pNH4])
    
        
    def Solo_M3pNa(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                    valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                    valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH])
    
    def Solo_M3pK(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                    valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                    valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH])
        
    def Solo_M4pK(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                    valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                    valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                    valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                    valid_M3pKpCH3CN, valid_M3pKpCH3OH])
    
    def Solo_M4pH(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                    valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                    valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                    valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                    valid_M3pKpCH3CN, valid_M3pKpCH3OH])
    
    def Solo_M4pNH4(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                    valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                    valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                    valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                    valid_M3pKpCH3CN, valid_M3pKpCH3OH])
    
    def Solo_M4pNa(ion_idx, mgf_file):
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
    
        valid_H = peaks.between(mz_H - prec_mass_error, mz_H + prec_mass_error, inclusive = True).sum() > 0
        valid_NH4 = peaks.between(mz_NH4 - prec_mass_error, mz_NH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_Na = peaks.between(mz_Na - prec_mass_error, mz_Na + prec_mass_error, inclusive = True).sum() > 0
        valid_K = peaks.between(mz_K - prec_mass_error, mz_K + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3CN = peaks.between(mz_HpCH3CN - prec_mass_error, mz_HpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_HpCH3OH = peaks.between(mz_HpCH3OH - prec_mass_error, mz_HpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3CN = peaks.between(mz_NapCH3CN - prec_mass_error, mz_NapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_NapCH3OH = peaks.between(mz_NapCH3OH - prec_mass_error, mz_NapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3CN = peaks.between(mz_KpCH3CN - prec_mass_error, mz_KpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_KpCH3OH = peaks.between(mz_KpCH3OH - prec_mass_error, mz_KpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M2pH = peaks.between(mz_M2pH - prec_mass_error, mz_M2pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNH4 = peaks.between(mz_M2pNH4 - prec_mass_error, mz_M2pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNa = peaks.between(mz_M2pNa - prec_mass_error, mz_M2pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pK = peaks.between(mz_M2pK - prec_mass_error, mz_M2pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3CN = peaks.between(mz_M2pHpCH3CN - prec_mass_error, mz_M2pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pHpCH3OH = peaks.between(mz_M2pHpCH3OH - prec_mass_error, mz_M2pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3CN = peaks.between(mz_M2pNapCH3CN - prec_mass_error, mz_M2pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pNapCH3OH = peaks.between(mz_M2pNapCH3OH - prec_mass_error, mz_M2pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3CN = peaks.between(mz_M2pKpCH3CN - prec_mass_error, mz_M2pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M2pKpCH3OH = peaks.between(mz_M2pKpCH3OH - prec_mass_error, mz_M2pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
    
        valid_M3pH = peaks.between(mz_M3pH - prec_mass_error, mz_M3pH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNH4 = peaks.between(mz_M3pNH4 - prec_mass_error, mz_M3pNH4 + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNa = peaks.between(mz_M3pNa - prec_mass_error, mz_M3pNa + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pK = peaks.between(mz_M3pK - prec_mass_error, mz_M3pK + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3CN = peaks.between(mz_M3pHpCH3CN - prec_mass_error, mz_M3pHpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pHpCH3OH = peaks.between(mz_M3pHpCH3OH - prec_mass_error, mz_M3pHpCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3CN = peaks.between(mz_M3pNapCH3CN - prec_mass_error, mz_M3pNapCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pNapCH3OH = peaks.between(mz_M3pNapCH3OH - prec_mass_error, mz_M3pNapCH3OH + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3CN = peaks.between(mz_M3pKpCH3CN - prec_mass_error, mz_M3pKpCH3CN + prec_mass_error, inclusive = True).sum() > 0
        valid_M3pKpCH3OH = peaks.between(mz_M3pKpCH3OH - prec_mass_error, mz_M3pKpCH3OH + prec_mass_error, inclusive = True).sum() > 0
            
        return sum([valid_H, valid_NH4, valid_Na, valid_K, valid_HpCH3CN, valid_HpCH3OH, valid_NapCH3CN, valid_NapCH3OH, valid_KpCH3CN, valid_KpCH3OH,
                    valid_M2pH, valid_M2pNH4, valid_M2pNa, valid_M2pK, valid_M2pHpCH3CN, valid_M2pHpCH3OH, valid_M2pNapCH3CN,
                    valid_M2pNapCH3OH, valid_M2pKpCH3CN, valid_M2pKpCH3OH, valid_M3pH, valid_M3pNH4, valid_M3pNa,
                    valid_M3pK, valid_M3pHpCH3CN, valid_M3pHpCH3OH, valid_M3pNapCH3CN, valid_M3pNapCH3OH,
                    valid_M3pKpCH3CN, valid_M3pKpCH3OH])
        
    def Not_validable(ion_idx, mgf_file):
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

    def Column_correction(table):
        drop_col = [i for i in table.columns if "Unnamed" in i]
        table.drop(drop_col, axis = 1, inplace = True)
        return table
    
    def Neutral_table(ion1_mz, adduct_table):
        """Computes all possible molecular masses (hypothetical neutrals) for an ion
        (ion 1) given the different ion species available in the adducts table.
        """
        neutral_table = adduct_table.copy()
        neutral_table['neutral_mass'] = [0.0]*len(neutral_table)
        for i in neutral_table.index :
            neutral_table.loc[i, 'neutral_mass'] = Neutral_mass_calculator(ion1_mz, 
                             neutral_table['Adduct_mass'][i], 
                             neutral_table['Mol_multiplier'][i], 
                             neutral_table['Charge'][i])
        return neutral_table
    
    def Ion_table(neutral_mass, adduct_table):
        """Computes all possible ion m/z values for an given molecular mass
        given the different ion species available in the adducts table.
        """        
        ion_table = adduct_table.copy()
        mz_list = []
        for i in ion_table.index:
            mz_list.append(Neutral_to_adduct(neutral_mass,
                                             ion_table.loc[i, "Adduct_mass"],
                                             ion_table.loc[i, "Mol_multiplier"],
                                             ion_table.loc[i, "Charge"]))
        ion_table['ion_mz'] = mz_list
        return ion_table

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
    
    def Samplewise_export(neg_csv_file, pos_csv_file, out_path, merged_edge_table, merged_node_table) : 
        print("Exporting sample-wise tables...")
        neg_csv = pd.read_csv(neg_csv_file, index_col ="row ID")
        pos_csv = pd.read_csv(pos_csv_file, index_col ="row ID")
        neg_csv = Column_correction(neg_csv)
        pos_csv = Column_correction(pos_csv)
        neg_csv.columns = neg_csv.columns.str.replace(".mzXML Peak area", "").str.replace('NEG_', '')
        pos_csv.columns = pos_csv.columns.str.replace(".mzXML Peak area", "").str.replace('POS_', '')
        neg_csv.drop(["row m/z", "row retention time"], axis = 1, inplace = True)
        pos_csv.drop(["row m/z", "row retention time"], axis = 1, inplace = True)
        samples = list(set(list(neg_csv.columns) + list(pos_csv.columns)))
        samples.sort()
        
        for sample in tqdm(samples):
            #sample = samples[0]
            ion_ids_neg = neg_csv.index[neg_csv[sample] > 0.0]
            ion_ids_pos = pos_csv.index[pos_csv[sample] > 0.0]
            
            #convert feature_ids to the new indexes
            tmp_table = merged_node_table[merged_node_table['status'] != "neg_neutral"]
            tmp_table = tmp_table[tmp_table['status'] != "pos_neutral"]
            tmp_table = tmp_table[tmp_table['status'] != "mix_neutral"]
            tmp_table_pos = tmp_table[tmp_table['ion_mode'] == "POS"]
            tmp_table_neg = tmp_table[tmp_table['ion_mode'] == "NEG"]
            ion_idx_neg = pd.Series(tmp_table_neg.index, index = tmp_table_neg['feature_id'])
            ion_idx_neg = list(ion_idx_neg[ion_ids_neg])
            ion_idx_pos = pd.Series(tmp_table_pos.index, index = tmp_table_pos['feature_id'])
            ion_idx_pos = list(ion_idx_pos[ion_ids_pos])
            ion_idx_mix = ion_idx_neg + ion_idx_pos
        
            # Get sample neutrals
            neutral_edges = merged_edge_table.loc[merged_edge_table["Adnotation"].dropna().index]
            kept_edges = [i for i in neutral_edges.index if neutral_edges.loc[i, "node_2"] in ion_idx_mix]
    
        
            # Get ion edges
            ion_edges = merged_edge_table[merged_edge_table['status'] != "neg_add_edge"]
            ion_edges = ion_edges[ion_edges['status'] != "pos_add_edge"]
            for i in ion_edges.index:
                if ion_edges.loc[i, "node_1"] in ion_idx_mix:
                    if ion_edges.loc[i, "node_2"] in ion_idx_mix:
                        kept_edges.append(i)
            kept_edges.sort()
            sample_edges = merged_edge_table.loc[kept_edges]
            sample_edges.sort_values('node_1', inplace = True)
            sample_edges.reset_index(inplace = True, drop = True)
            
            kept_nodes = list(set(list(sample_edges['node_1']) + list(sample_edges['node_2'])))
            kept_nodes.sort()
            sample_nodes = merged_node_table.loc[kept_nodes].copy()
            sample_nodes.drop(pd.Series(samples) + ".mzXML Peak area", axis = 1, inplace = True)
            sample_nodes[sample] = merged_node_table[sample + ".mzXML Peak area"]
            
            sample_nodes.to_csv(out_path + "MIX_" + sample + "_nodes.csv", index_label = "Index")
            sample_edges.to_csv(out_path + "MIX_" + sample + "_edges.csv", index_label = "Index")
        return
    
    

    import os
    import pandas as pd
    from pandas.core.common import flatten
    import sys
    from tqdm import tqdm
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    def Spectrum_processing(s):
        s = default_filters(s)
        return s

    # Load parameters:
    in_path_full_neg= params['neg_out_0']
    in_path_full_pos= params['pos_out_0']
    neg_mgf= params['neg_mgf']
    pos_mgf= params['pos_mgf']
    pos_csv= params['pos_csv']
    neg_csv= params['neg_csv']
    adnotation_pos_full= params['pos_out_3_1']
    adnotation_neg_full= params['neg_out_3_1']
    out_full= params['mix_out_4_1']
    out_samples= params['mix_out_4_2']
    neg_nodes_full = "NEG_full_nodes.csv"
    neg_edges_full = "NEG_full_edges.csv"
    pos_nodes_full = "POS_full_nodes.csv"
    pos_edges_full = "POS_full_edges.csv"
    
    adduct_table_primary_neg= params['mm_addtable_primary_neg']
    adduct_table_secondary_neg= params['mm_addtable_secondary_neg']
    adduct_table_primary_pos= params['mm_addtable_primary_pos']
    adduct_table_secondary_pos= params['mm_addtable_secondary_pos']
    mass_error= params['mm_mass_error']
    prec_mass_error = params['mm_prec_mass_error']
    rt_error= params['mm_rt_error']  
    bnr_list_neg = params['mm_bnr_neg']
    bnr_list_pos = params['mm_bnr_pos']
    
    #Load and filter MGF files
    print("Loading and filtering NEG MGF file...")
    mgf_file_neg = list(load_from_mgf(in_path_full_neg + neg_mgf))
    mgf_file_neg = [Spectrum_processing(s) for s in mgf_file_neg]
    print("Loading and filtering POS MGF file...")
    mgf_file_pos = list(load_from_mgf(in_path_full_pos + pos_mgf))
    mgf_file_pos = [Spectrum_processing(s) for s in mgf_file_pos]
    
    
    #Load adduct tables
    adduct_table_primary_neg = pd.read_csv("./params/" + adduct_table_primary_neg, sep = "\t")
    adduct_table_secondary_neg = pd.read_csv("./params/" + adduct_table_secondary_neg, sep = "\t")
    adduct_table_merged_neg = adduct_table_primary_neg.append(adduct_table_secondary_neg, ignore_index = True)
    adduct_table_primary_pos = pd.read_csv("./params/" + adduct_table_primary_pos, sep = "\t")
    adduct_table_secondary_pos = pd.read_csv("./params/" + adduct_table_secondary_pos, sep = "\t")
    adduct_table_merged_pos = adduct_table_primary_pos.append(adduct_table_secondary_pos, ignore_index = True)
    
    # Produce base neutral adduct tables
    adduct_table_base_neg = [adduct_table_merged_neg.index[adduct_table_merged_neg['Adduct'] == adduct][0] for adduct in bnr_list_neg]
    adduct_table_base_neg = adduct_table_merged_neg.loc[adduct_table_base_neg]
    adduct_table_base_pos = [adduct_table_merged_pos.index[adduct_table_merged_pos['Adduct'] == adduct][0] for adduct in bnr_list_pos]
    adduct_table_base_pos = adduct_table_merged_pos.loc[adduct_table_base_pos]
    
    # Loading files: POS and NEG node tables, edge tables and MZmine CSV output
    neg_mzmine_csv = pd.read_csv(in_path_full_neg + neg_csv, index_col = "row ID")
    neg_mzmine_csv = Column_correction(neg_mzmine_csv)
    pos_mzmine_csv = pd.read_csv(in_path_full_pos + pos_csv, index_col = "row ID")
    pos_mzmine_csv = Column_correction(pos_mzmine_csv)
    edge_table_all_neg = pd.read_csv(adnotation_neg_full + neg_edges_full, index_col = "Index")
    node_table_all_neg = pd.read_csv(adnotation_neg_full + neg_nodes_full, index_col = "feature_id")
    edge_table_all_pos = pd.read_csv(adnotation_pos_full + pos_edges_full, index_col = "Index")
    node_table_all_pos = pd.read_csv(adnotation_pos_full + pos_nodes_full, index_col = "feature_id")
    
    # Setting ion modes in tables
    node_table_all_neg['ion_mode'] = ['NEG']*len(node_table_all_neg)
    edge_table_all_neg['ion_mode'] = ['NEG']*len(edge_table_all_neg)
    node_table_all_pos['ion_mode'] = ['POS']*len(node_table_all_pos)
    edge_table_all_pos['ion_mode'] = ['POS']*len(edge_table_all_pos)
    
    # New indexes will be used when merging POS and NEG tables, feature_IDs are saved in a dedicated column
    node_table_all_neg['feature_id'] = node_table_all_neg.index
    node_table_all_pos['feature_id'] = node_table_all_pos.index
    
    # Set samples in node tables :
    neg_mzmine_csv.drop(['row m/z', 'row retention time'], axis = 1, inplace = True)
    pos_mzmine_csv.drop(['row m/z', 'row retention time'], axis = 1, inplace = True)
    neg_mzmine_csv.columns = pd.Series(neg_mzmine_csv.columns).str.replace('NEG_', '')
    pos_mzmine_csv.columns = pd.Series(pos_mzmine_csv.columns).str.replace('POS_', '')
    node_table_all_neg["samples"] = ['']*len(node_table_all_neg)
    node_table_all_pos["samples"] = ['']*len(node_table_all_pos)
    neg_ions_idx = node_table_all_neg[node_table_all_neg['status'] != "neutral"].index
    pos_ions_idx = node_table_all_pos[node_table_all_pos['status'] != "neutral"].index
    print("Adding samples to neg ions...")
    for ion in tqdm(neg_ions_idx):
        samples = list(neg_mzmine_csv.columns[neg_mzmine_csv.loc[ion] > 0])
        node_table_all_neg.loc[ion, "samples"] = '|'.join(samples)
    print("Adding samples to pos ions...")
    for ion in tqdm(pos_ions_idx):
        samples = list(pos_mzmine_csv.columns[pos_mzmine_csv.loc[ion] > 0])
        node_table_all_pos.loc[ion, "samples"] = '|'.join(samples)
    
    # Add samples to neutrals
    neg_neutrals_idx = node_table_all_neg[node_table_all_neg['status'] == "neutral"].index
    pos_neutrals_idx = node_table_all_pos[node_table_all_pos['status'] == "neutral"].index
    print("Adding samples to neg neutrals...")
    for neutral in tqdm(neg_neutrals_idx):
        tmp_ions = list(edge_table_all_neg['node_2'][edge_table_all_neg['node_1'] == neutral])
        tmp_samples = []
        for ion in tmp_ions:
            tmp_samples += node_table_all_neg.loc[ion, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        node_table_all_neg.loc[neutral, "samples"] = tmp_samples
    print("Adding samples to pos neutrals...")
    for neutral in tqdm(pos_neutrals_idx):
        tmp_ions = list(edge_table_all_pos['node_2'][edge_table_all_pos['node_1'] == neutral])
        tmp_samples = []
        for ion in tmp_ions:
            tmp_samples += node_table_all_pos.loc[ion, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        node_table_all_pos.loc[neutral, "samples"] = tmp_samples
    
    
    # Produce the merged mode edge and node tables
    node_table_ions_neg= node_table_all_neg[node_table_all_neg['status'] != "neutral"]
    node_table_ions_pos= node_table_all_pos[node_table_all_pos['status'] != "neutral"]
    merged_node_table = pd.DataFrame()
    merged_node_table = merged_node_table.append(node_table_ions_neg, ignore_index = True)
    merged_node_table = merged_node_table.append(node_table_ions_pos, ignore_index = True)
    merged_edge_table = pd.DataFrame()
    merged_edge_table = merged_edge_table.append(edge_table_all_neg, ignore_index = True)
    merged_edge_table = merged_edge_table.append(edge_table_all_pos, ignore_index = True)
    
    #Chekc from NEG to POS if there are neutral nodes that might match
    neg_neutrals = node_table_all_neg[node_table_all_neg['status'] == "neutral"].copy()
    neg_neutrals['feature_id'] = [None]*len(neg_neutrals)
    pos_neutrals = node_table_all_pos[node_table_all_pos['status'] == "neutral"].copy()
    pos_neutrals['feature_id'] = [None]*len(pos_neutrals)
    merged_neutrals = []
    for neutral in tqdm(neg_neutrals_idx): # NEG -> POS
        #neutral = neg_neutrals_idx[0]
        mz_neg = neg_neutrals.loc[neutral, "mz"]
        rt_neg = neg_neutrals.loc[neutral, "rt"]
        samples_neg = set(neg_neutrals.loc[neutral, "samples"].split('|'))
        counter_pos = pos_neutrals[pos_neutrals["mz"].between(mz_neg - mass_error, mz_neg + mass_error, inclusive = True)]
        counter_pos = counter_pos[counter_pos["rt"].between(rt_neg - rt_error, rt_neg + rt_error, inclusive = True)] 
        if len(counter_pos) == 0 : continue
        counter_pos['d_rt'] = abs(counter_pos['rt']-rt_neg)
        shared_samples = []
        for i in counter_pos.index:
            samples_pos = counter_pos.loc[i, "samples"].split('|')
            shared_samples.append(len(samples_neg.intersection(samples_pos)))
        counter_pos['shared_samples'] = shared_samples
        counter_pos = counter_pos[counter_pos['shared_samples'] > 0]
        counter_pos = counter_pos[counter_pos['shared_samples'] == counter_pos['shared_samples'].max()]
        counter_pos = counter_pos[counter_pos['d_rt'] == counter_pos['d_rt'].min()]
        if len(counter_pos) > 1 : sys.exit("MORE THAN 1 COUTNER NEUTRAL AT SAME RT")
        elif len(counter_pos) == 0 : continue
        else:
            mz_pos = counter_pos["mz"].iloc[0]
            rt_pos = counter_pos["rt"].iloc[0]
            merged_neutrals.append((neutral, mz_neg, rt_neg, counter_pos.index[0], mz_pos, rt_pos))
    
    merged_neutrals = pd.DataFrame(merged_neutrals, columns = ['neg_neutral', 'mass_neg',
                                                               'rt_neg', 'pos_neutral',
                                                               'mass_pos', 'rt_pos'])
    merged_neutrals['d_rt'] = abs(merged_neutrals['rt_pos'] - merged_neutrals['rt_neg'])
    
    # Sometimes, several NEG neutrals might match the same POS neutral. 
    # The NEG neutral with the closest RT to the POS neutral is chosen
    for neutral in merged_neutrals['pos_neutral'].unique():
        tmp_table = merged_neutrals[merged_neutrals['pos_neutral'] == neutral].copy()
        if len(tmp_table) > 1 : 
            tmp_table.sort_values('d_rt', inplace = True)
            for i in tmp_table.index[1:]:
                merged_neutrals.drop(i, inplace = True)
    
    # A transitions table is produced containing the neutrals to be merged and their 
    # new ID. Later will also contain non-merged neutrals and their new IDs
    transitions_table = list()
    print("Adding merged neutrals to the merged node table...")
    for i in tqdm(merged_neutrals.index):
        new_idx = merged_node_table.index.max() + 1 
        neutral_neg = merged_neutrals.loc[i, "neg_neutral"]
        neutral_pos = merged_neutrals.loc[i, "pos_neutral"]
        neg_ions = list(edge_table_all_neg['node_2'][edge_table_all_neg['node_1'] == neutral_neg])
        neg_count = len(neg_ions)
        neg_mz = node_table_all_neg.loc[neutral_neg, "mz"]
        neg_rt = node_table_all_neg.loc[neutral_neg, "rt"]
        neg_tic = node_table_all_neg.loc[neutral_neg, "TIC"]
        neg_samples = node_table_all_neg.loc[neutral_neg, "samples"].split('|')
        pos_ions = list(edge_table_all_pos['node_2'][edge_table_all_pos['node_1'] == neutral_pos])
        pos_count = len(pos_ions)
        pos_mz = node_table_all_pos.loc[neutral_pos, "mz"]
        pos_rt = node_table_all_pos.loc[neutral_pos, "rt"]
        pos_tic = node_table_all_pos.loc[neutral_pos, "TIC"]
        pos_samples = node_table_all_pos.loc[neutral_pos, "samples"].split('|')
        mix_mz = round(((neg_mz*neg_count) + (pos_mz*pos_count))/(neg_count + pos_count), 4)
        mix_rt = round(((neg_rt*neg_count) + (pos_rt*pos_count))/(neg_count + pos_count), 3)
        mix_tic = neg_tic + pos_tic
        mix_samples = list(set(neg_samples + pos_samples))
        mix_samples.sort()
        mix_samples = '|'.join(mix_samples)
        merged_node_table.loc[new_idx] = [mix_mz, mix_rt, mix_tic, 0, None, 0, "neutral",
                             len(neg_ions) + len(pos_ions), None, "MIX", None, mix_samples]
        transitions_table.append((new_idx, neutral_neg, neutral_pos))
    transitions_table = pd.DataFrame(transitions_table, columns = ['MIX', 'NEG', 'POS'])   
    
    # Drop neutrals that were merged and redo the process for the remaining neutrals
    # that might have escaped, this time in POS to NEG for a change.
    neg_neutrals.drop(transitions_table['NEG'], inplace = True)
    pos_neutrals.drop(transitions_table['POS'], inplace = True)
    merged_neutrals = []
    for neutral in tqdm(pos_neutrals.index): # POS -> NEG
        #neutral = pos_neutrals.index[0]
        mz_pos = pos_neutrals.loc[neutral, "mz"]
        rt_pos = pos_neutrals.loc[neutral, "rt"]
        samples_pos = set(pos_neutrals.loc[neutral, "samples"].split('|'))
        counter_neg = neg_neutrals[neg_neutrals["mz"].between(mz_pos - mass_error, mz_pos + mass_error, inclusive = True)]
        counter_neg = counter_neg[counter_neg["rt"].between(rt_pos - rt_error, rt_pos + rt_error, inclusive = True)] 
        if len(counter_neg) == 0 : continue
        counter_neg['d_rt'] = abs(counter_neg['rt']-rt_pos)
        shared_samples = []
        for i in counter_neg.index:
            samples_neg = counter_neg.loc[i, "samples"].split('|')
            shared_samples.append(len(samples_pos.intersection(samples_neg)))
        counter_neg['shared_samples'] = shared_samples
        counter_neg = counter_neg[counter_neg['shared_samples'] > 0]
        counter_neg = counter_neg[counter_neg['shared_samples'] == counter_neg['shared_samples'].max()]
        counter_neg = counter_neg[counter_neg['d_rt'] == counter_neg['d_rt'].min()]
        if len(counter_neg) > 1 : sys.exit("MORE THAN 1 COUTNER NEUTRAL AT SAME RT")
        elif len(counter_neg) == 0 : continue
        else:
            mz_neg = counter_neg["mz"].iloc[0]
            rt_neg = counter_neg["rt"].iloc[0]
            merged_neutrals.append((counter_neg.index[0], mz_neg, rt_neg, neutral, mz_pos, rt_pos))
    merged_neutrals = pd.DataFrame(merged_neutrals, columns = ['neg_neutral', 'mass_neg',
                                                               'rt_neg', 'pos_neutral',
                                                               'mass_pos', 'rt_pos'])
    merged_neutrals['d_rt'] = abs(merged_neutrals['rt_pos'] - merged_neutrals['rt_neg'])
    
    # Same as before, if there are several POS neutrals per NEG neutral, choose the
    # closest one by RT
    for neutral in merged_neutrals['pos_neutral'].unique():
        tmp_table = merged_neutrals[merged_neutrals['pos_neutral'] == neutral].copy()
        if len(tmp_table) > 1 : 
            tmp_table.sort_values('d_rt', inplace = True)
            for i in tmp_table.index[1:]:
                merged_neutrals.drop(i, inplace = True)
    
    # Add the new merged neutrals to the transitions table
    print("Adding merged neutrals to the merged node table...")
    for i in tqdm(merged_neutrals.index):
        new_idx = merged_node_table.index.max() + 1 
        neutral_neg = merged_neutrals.loc[i, "neg_neutral"]
        neutral_pos = merged_neutrals.loc[i, "pos_neutral"]
        neg_ions = list(edge_table_all_neg['node_2'][edge_table_all_neg['node_1'] == neutral_neg])
        neg_count = len(neg_ions)
        neg_mz = node_table_all_neg.loc[neutral_neg, "mz"]
        neg_rt = node_table_all_neg.loc[neutral_neg, "rt"]
        neg_tic = node_table_all_neg.loc[neutral_neg, "TIC"]
        neg_samples = node_table_all_neg.loc[neutral_neg, "samples"].split('|')
        pos_ions = list(edge_table_all_pos['node_2'][edge_table_all_pos['node_1'] == neutral_pos])
        pos_count = len(pos_ions)
        pos_mz = node_table_all_pos.loc[neutral_pos, "mz"]
        pos_rt = node_table_all_pos.loc[neutral_pos, "rt"]
        pos_tic = node_table_all_pos.loc[neutral_pos, "TIC"]
        pos_samples = node_table_all_pos.loc[neutral_pos, "samples"].split('|')
        mix_mz = round(((neg_mz*neg_count) + (pos_mz*pos_count))/(neg_count + pos_count), 4)
        mix_rt = round(((neg_rt*neg_count) + (pos_rt*pos_count))/(neg_count + pos_count), 3)
        mix_tic = neg_tic + pos_tic
        mix_samples = list(set(neg_samples + pos_samples))
        mix_samples.sort()
        mix_samples = '|'.join(mix_samples)
        merged_node_table.loc[new_idx] = [mix_mz, mix_rt, mix_tic, 0, None, 0, "neutral",
                             len(neg_ions) + len(pos_ions), None, "MIX", None, mix_samples]
        transitions_table.loc[transitions_table.index.max() + 1] = [new_idx, neutral_neg, neutral_pos]

    # Add non merged neutrals to the merged_node_table (using their new IDs) and 
    # to the transitions table
    # First, eliminate merged_neutrals 
    intersect_neutrals = set(neg_neutrals.index)
    intersect_neutrals = intersect_neutrals.intersection(transitions_table['NEG'])
    neg_neutrals.drop(intersect_neutrals, inplace = True)
    
    intersect_neutrals = set(pos_neutrals.index)
    intersect_neutrals = intersect_neutrals.intersection(transitions_table['POS'])    
    pos_neutrals.drop(intersect_neutrals, inplace = True)
    
    print('Adding non-merged neg neutrals to the merged node table...')
    for neutral in tqdm(neg_neutrals.index):
        new_idx = merged_node_table.index.max() + 1 
        merged_node_table.loc[new_idx] = neg_neutrals.loc[neutral].copy()
        new_row = pd.Series([new_idx, neutral, None], index = ["MIX", "NEG", "POS"])
        transitions_table = transitions_table.append(new_row, ignore_index = True)

    
    print('Adding non-merged pos neutrals to the merged node table...')
    for neutral in tqdm(pos_neutrals.index):
        new_idx = merged_node_table.index.max() + 1 
        merged_node_table.loc[new_idx] = pos_neutrals.loc[neutral].copy()
        new_row = pd.Series([new_idx, None, neutral], index = ["MIX", "NEG", "POS"])
        transitions_table = transitions_table.append(new_row, ignore_index = True)
    
    # Replace old IDs by new ones in the merged edge table.
    print('Resetting edge table IDs...')
    for i in tqdm(merged_edge_table.index):
        ion_mode = merged_edge_table.loc[i, "ion_mode"]
        if ion_mode == "NEG" : opposed_mode = "POS"
        else : opposed_mode = "NEG"
        merged_node_table_sub = merged_node_table[merged_node_table['ion_mode'] != opposed_mode]
        if merged_edge_table.loc[i, "status"] == "add_edge":
            node_1 = merged_edge_table.loc[i, "node_1"]
            node_2 = merged_edge_table.loc[i, "node_2"]
            new_node_1 = int(transitions_table["MIX"][transitions_table[ion_mode] == node_1].iloc[0])
            new_node_2 = merged_node_table_sub.index[merged_node_table_sub["feature_id"] == node_2][0]
            merged_edge_table.loc[i, "node_1"] = new_node_1
            merged_edge_table.loc[i, "node_2"] = new_node_2
        else:
            node_1 = merged_edge_table.loc[i, "node_1"]
            node_2 = merged_edge_table.loc[i, "node_2"]
            new_node_1 = merged_node_table_sub.index[merged_node_table_sub["feature_id"] == node_1][0]
            new_node_2 = merged_node_table_sub.index[merged_node_table_sub["feature_id"] == node_2][0]
            merged_edge_table.loc[i, "node_1"] = new_node_1
            merged_edge_table.loc[i, "node_2"] = new_node_2
    
    # Search for remaining ion annotations
    remains_table_pos = merged_node_table[merged_node_table['status'] != "neutral"]
    remains_table_pos = remains_table_pos[remains_table_pos['status'] != "adduct"] 
    remains_table_pos = remains_table_pos[remains_table_pos['ion_mode'] == "POS"]
    remains_table_neg = merged_node_table[merged_node_table['status'] != "neutral"]
    remains_table_neg = remains_table_neg[remains_table_neg['status'] != "adduct"] 
    remains_table_neg = remains_table_neg[remains_table_neg['ion_mode'] == "NEG"]
    node_table_neutrals_pos = merged_node_table[merged_node_table["status"] == "neutral"]
    node_table_neutrals_pos = node_table_neutrals_pos[node_table_neutrals_pos["ion_mode"] == "POS"]
    node_table_neutrals_neg = merged_node_table[merged_node_table["status"] == "neutral"]
    node_table_neutrals_neg = node_table_neutrals_neg[node_table_neutrals_neg["ion_mode"] == "NEG"]      

    print('Linking NEG neutrals to single POS ions...')
    candidates_table = list()
    for i in tqdm(node_table_neutrals_neg.index):
        mol_mass = node_table_neutrals_neg.loc[i, "mz"]
        mol_rt = node_table_neutrals_neg.loc[i, "rt"]
        mol_samples = set(node_table_neutrals_neg.loc[i, "samples"].split('|'))
        ion_table = Ion_table(mol_mass, adduct_table_merged_pos)
        ion_hits = list()
        hit_table_rt = remains_table_pos[remains_table_pos['rt'].between(mol_rt - rt_error, mol_rt + rt_error, inclusive = True)].copy()
        
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(mol_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        
        for j in ion_table.index:
            ion_mz = ion_table.loc[j, "ion_mz"]
            hit_table = hit_table_rt[hit_table_rt['mz'].between(ion_mz - mass_error, ion_mz + mass_error, inclusive = True)].copy()
            if len(hit_table) >0 :
                ion_hits.append('|'.join(hit_table.index.astype(str)))
            else:
                ion_hits.append(None)
        ion_table['ion_hits'] = ion_hits
        ion_table = ion_table[~ion_table['ion_hits'].isnull()]
        for j in ion_table.index:
            adduct = ion_table.loc[j, "Adduct"]
            for k in ion_table.loc[j, "ion_hits"].split('|'):
                k = int(k)
                d_rt = abs(mol_rt - remains_table_pos.loc[k, 'rt'])
                d_mz = abs(mol_mass - remains_table_pos.loc[k, 'mz'])
                candidates_table.append((k, i, adduct, d_rt, d_mz))
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Migration of ions:
    unique_hits = candidates_table['ion_idx'].unique().tolist()
    selected_neutrals = list()
    selected_adducts = list()
    delta_rts = list()
    delta_mzs = list()
    for i in unique_hits:
        tmp_table = candidates_table[candidates_table['ion_idx'] == i]
        selected_neutral = tmp_table['d_rt'].idxmin()
        selected_neutrals.append(tmp_table.loc[selected_neutral, "neutral"])
        selected_adducts.append(tmp_table.loc[selected_neutral, "adduct"])
        delta_rts.append(tmp_table.loc[selected_neutral, "d_rt"])
        delta_mzs.append(tmp_table.loc[selected_neutral, "d_mz"])
    candidates_table= list(zip(unique_hits, selected_neutrals, selected_adducts, delta_rts, delta_mzs))
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Report results and node and edge tables:
    print('Updating node and edge tables...')
    for i in tqdm(candidates_table.index):
        # Retrieve data for ion and neutral
        ion_idx = candidates_table.loc[i, "ion_idx"]
        neutral_idx = candidates_table.loc[i, "neutral"]
        adduct = candidates_table.loc[i, "adduct"]
        d_rt = candidates_table.loc[i, "d_rt"]
        d_mz = candidates_table.loc[i, "d_mz"]
        mgf_idx = int(merged_node_table.loc[ion_idx, "mgf_index"])
        adduct_code = adduct_table_merged_pos['Adduct_code'][adduct_table_merged_pos['Adduct'] == adduct].iloc[0]
        Species_rule = Validator_choice(adduct_code, "POS")
        rule_points = Species_rule(mgf_idx, mgf_file_pos)
        ion_tic = merged_node_table.loc[ion_idx, 'TIC']
        
        # Update node table with ion data
        merged_node_table.loc[ion_idx, 'status'] = "adduct"
        merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
        merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
        
        # Update node table with neutral data
        merged_node_table.loc[neutral_idx, "TIC"] += ion_tic
        merged_node_table.loc[neutral_idx, "adduct_count"] += 1
        merged_node_table.loc[neutral_idx, "ion_mode"] = "MIX"
        tmp_samples = merged_node_table.loc[neutral_idx, "samples"].split('|')
        tmp_samples += merged_node_table.loc[ion_idx, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        merged_node_table.loc[neutral_idx, "samples"] = tmp_samples

        # Update edge table :
        del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
        del_edge = del_edge.index[del_edge['node_2'] == ion_idx]
        if len(del_edge) > 0 : merged_edge_table.drop(del_edge[0], inplace = True)
        new_idx = merged_edge_table.index.max() + 1 
        merged_edge_table.loc[new_idx] = [neutral_idx, ion_idx, 0, 0, 0, d_rt,
                             d_mz, "add_edge", None, adduct, adduct, "POS"]
    
    print('Linking POS neutrals to single NEG ions...')
    candidates_table = list()
    for i in tqdm(node_table_neutrals_pos.index):
        mol_mass = node_table_neutrals_pos.loc[i, "mz"]
        mol_rt = node_table_neutrals_pos.loc[i, "rt"]
        mol_samples = set(node_table_neutrals_pos.loc[i, "samples"].split('|'))
        ion_table = Ion_table(mol_mass, adduct_table_merged_neg)
        ion_hits = list()
        hit_table_rt = remains_table_neg[remains_table_neg['rt'].between(mol_rt - rt_error, mol_rt + rt_error, inclusive = True)].copy()
        
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(mol_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        
        for j in ion_table.index:
            ion_mz = ion_table.loc[j, "ion_mz"]
            hit_table = hit_table_rt[hit_table_rt['mz'].between(ion_mz - mass_error, ion_mz + mass_error, inclusive = True)].copy()
            if len(hit_table) >0 :
                ion_hits.append('|'.join(hit_table.index.astype(str)))
            else:
                ion_hits.append(None)
        ion_table['ion_hits'] = ion_hits
        ion_table = ion_table[~ion_table['ion_hits'].isnull()]
        for j in ion_table.index:
            adduct = ion_table.loc[j, "Adduct"]
            for k in ion_table.loc[j, "ion_hits"].split('|'):
                k = int(k)
                d_rt = abs(mol_rt - remains_table_neg.loc[k, 'rt'])
                d_mz = abs(mol_mass - remains_table_neg.loc[k, 'mz'])
                candidates_table.append((k, i, adduct, d_rt, d_mz))
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Migration of ions:
    unique_hits = candidates_table['ion_idx'].unique().tolist()
    selected_neutrals = list()
    selected_adducts = list()
    delta_rts = list()
    delta_mzs = list()
    for i in unique_hits:
        tmp_table = candidates_table[candidates_table['ion_idx'] == i]
        selected_neutral = tmp_table['d_rt'].idxmin()
        selected_neutrals.append(tmp_table.loc[selected_neutral, "neutral"])
        selected_adducts.append(tmp_table.loc[selected_neutral, "adduct"])
        delta_rts.append(tmp_table.loc[selected_neutral, "d_rt"])
        delta_mzs.append(tmp_table.loc[selected_neutral, "d_mz"])
    candidates_table= list(zip(unique_hits, selected_neutrals, selected_adducts, delta_rts, delta_mzs))
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Report results and node and edge tables:
    print('Updating node and edge tables...')
    for i in tqdm(candidates_table.index):
        # Retrieve data for ion and neutral
        ion_idx = candidates_table.loc[i, "ion_idx"]
        neutral_idx = candidates_table.loc[i, "neutral"]
        adduct = candidates_table.loc[i, "adduct"]
        d_rt = candidates_table.loc[i, "d_rt"]
        d_mz = candidates_table.loc[i, "d_mz"]
        mgf_idx = int(merged_node_table.loc[ion_idx, "mgf_index"])
        adduct_code = adduct_table_merged_neg['Adduct_code'][adduct_table_merged_neg['Adduct'] == adduct].iloc[0]
        Species_rule = Validator_choice(adduct_code, "NEG")
        rule_points = Species_rule(mgf_idx, mgf_file_neg)
        ion_tic = merged_node_table.loc[ion_idx, 'TIC']
        
        # Update node table with ion data
        merged_node_table.loc[ion_idx, 'status'] = "adduct"
        merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
        merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
        
        # Update node table with neutral data
        merged_node_table.loc[neutral_idx, "TIC"] += ion_tic
        merged_node_table.loc[neutral_idx, "adduct_count"] += 1
        merged_node_table.loc[neutral_idx, "ion_mode"] = "MIX"
        tmp_samples = merged_node_table.loc[neutral_idx, "samples"].split('|')
        tmp_samples += merged_node_table.loc[ion_idx, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        merged_node_table.loc[neutral_idx, "samples"] = tmp_samples

        # Update edge table :
        del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
        del_edge = del_edge.index[del_edge['node_2'] == ion_idx]
        if len(del_edge) > 0 : merged_edge_table.drop(del_edge[0], inplace = True)
        new_idx = merged_edge_table.index.max() + 1 
        merged_edge_table.loc[new_idx] = [neutral_idx, ion_idx, 0, 0, 0, d_rt,
                             d_mz, "add_edge", None, adduct, adduct, "NEG"]
    
    # Produce cluster IDs. Making these is essential before opposed mode singletons
    # search, to filter out singletons already present in molecular clusters.
    node_pool = list(merged_node_table.index)
    singletons = list(merged_edge_table["node_1"][merged_edge_table['status'] == "self_edge"])
    node_pool = list(set(node_pool) - set(singletons))
    cluster_list = []
    cluster_size_list = []
    total_nodes = len(node_pool)
    while len(node_pool) > 0:
        new_cluster = [node_pool[0]]
        cluster_size = 0
        perc = round((1-(len(node_pool)/total_nodes))*100,1)
        sys.stdout.write("\rDefining new clusters : {0}%".format(perc))
        sys.stdout.flush()
        while cluster_size != len(new_cluster):
            cluster_size = len(new_cluster)
            tmp_idx = []
            for i in new_cluster:
                tmp_idx += list(merged_edge_table.index[merged_edge_table['node_1'] == i])
                tmp_idx += list(merged_edge_table.index[merged_edge_table['node_2'] == i])
            new_cluster += list(merged_edge_table.loc[tmp_idx, 'node_1'])
            new_cluster += list(merged_edge_table.loc[tmp_idx, 'node_2'])
            new_cluster = list(set(new_cluster))
        new_cluster.sort()
        node_pool = list(set(node_pool) - set(new_cluster))
        cluster_size_list.append(len(new_cluster))
        cluster_list.append('|'.join(list(map(str, new_cluster))))

    cluster_table= pd.DataFrame()
    cluster_table['cluster'] = cluster_list
    cluster_table['cluster_size'] = cluster_size_list
    cluster_table.sort_values('cluster_size', ascending = False, inplace = True)
    cluster_table.reset_index(drop = True, inplace = True)
    
    # Identify molecular clusters
    cluster_molecular = list()
    for i in cluster_table.index:
        node_list = cluster_table.loc[i, "cluster"].split('|')
        node_list = list(map(int, node_list))
        tmp_table_1 = merged_node_table.loc[node_list]
        if sum(tmp_table_1['status'] == "neutral") > 0 :
            cluster_molecular.append(True)
        else:
            cluster_molecular.append(False)
    cluster_table["molecular_cluster"] = cluster_molecular
    
    merged_node_table['cluster_id'] = [-1]*len(merged_node_table)
    print('Assigning new cluster indexes...')
    for i in tqdm(cluster_table.index):
        node_list = list(map(int, cluster_table.loc[i, 'cluster'].split('|')))
        for j in node_list :
            merged_node_table.loc[j, 'cluster_id'] = i
    
    # Connect singletons to singletons (only precursors and nodes from non molecular clusters)
    remains_table = cluster_table.index[~cluster_table["molecular_cluster"]].tolist()
    remains_table = [merged_node_table.index[merged_node_table['cluster_id'] == i].tolist() for i in remains_table]
    remains_table = list(flatten(remains_table)) # Added precursors and fragments from non-molecular clusters
    remains_table += merged_node_table.index[merged_node_table['cluster_id'] == -1].tolist() # Added singleton nodes
    precursor_ions = list(set(merged_node_table.index) - set(remains_table))
    precursor_ions = merged_node_table.loc[precursor_ions]
    precursor_ions = precursor_ions.index[precursor_ions['status'] == "precursor"].tolist()
    remains_table += precursor_ions # Added precursors from molecular clusters
    remains_table.sort()
    remains_table = merged_node_table.loc[remains_table]
    remains_table_pos = remains_table[remains_table['ion_mode'] == "POS"]
    remains_table_neg = remains_table[remains_table['ion_mode'] == "NEG"]

    print('Linking NEG singletons to POS singletons...')
    candidates_table = list()
    for i in tqdm(remains_table_neg.index):
        ion_1_mz = remains_table_neg.loc[i, "mz"]
        ion_1_rt = remains_table_neg.loc[i, "rt"]
        ion_1_samples = set(remains_table_neg.loc[i, "samples"].split('|'))
        neutral_table = Neutral_table(ion_1_mz, adduct_table_base_neg)
        hit_table_rt = remains_table_pos[remains_table_pos['rt'].between(ion_1_rt - rt_error, ion_1_rt + rt_error, inclusive = True)].copy()
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(ion_1_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        if len(hit_table_rt) == 0 : continue
        for j in neutral_table.index:
            ion_table = Ion_table(neutral_table.loc[j, "neutral_mass"], adduct_table_base_pos)
            for k in ion_table.index:
                hit_table = hit_table_rt[hit_table_rt['mz'].between(ion_table.loc[k, "ion_mz"] - mass_error, ion_table.loc[k, "ion_mz"] + mass_error, inclusive = True)]
                for l in hit_table.index:
                    candidates_table.append((i,
                                             neutral_table.loc[j, "Adduct"],
                                             neutral_table.loc[j, "Complexity"],
                                             l,
                                             ion_table.loc[k, "Adduct"],
                                             ion_table.loc[k, "Complexity"]))
    candidates_table = pd.DataFrame(candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity'])

    # Get the best hypotheses
    unique_negs = candidates_table['neg_ion'].unique().tolist()
    neg_ion_list = list()
    neg_adduct_list = list()
    neg_complexity_list = list()
    pos_ion_list = list()
    pos_adduct_list = list()
    pos_complexity_list = list()
    print('Selecting best NEG annotations...')
    for i in tqdm(unique_negs) :
        tmp_table_1 = candidates_table[candidates_table['neg_ion'] == i]
        unique_adducts = tmp_table_1['neg_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['neg_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['neg_complexity'].sum() + tmp_table_2['pos_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['neg_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            neg_ion_list.append(tmp_table_1.loc[j, "neg_ion"])
            neg_adduct_list.append(tmp_table_1.loc[j, "neg_adduct"])
            neg_complexity_list.append(tmp_table_1.loc[j, "neg_complexity"])
            pos_ion_list.append(tmp_table_1.loc[j, "pos_ion"])
            pos_adduct_list.append(tmp_table_1.loc[j, "pos_adduct"])
            pos_complexity_list.append(tmp_table_1.loc[j, "pos_complexity"])
    candidates_table = list(zip(neg_ion_list, neg_adduct_list, neg_complexity_list, pos_ion_list, pos_adduct_list, pos_complexity_list))
    candidates_table = pd.DataFrame(candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity'])

    # Resolve POS ions with multiple annotations
    unique_pos = candidates_table['pos_ion'].unique().tolist()
    neg_ion_list = list()
    neg_adduct_list = list()
    neg_complexity_list = list()
    pos_ion_list = list()
    pos_adduct_list = list()
    pos_complexity_list = list()
    print('Selecting best POS annotations...')
    for i in tqdm(unique_pos):
        tmp_table_1 = candidates_table[candidates_table['pos_ion'] == i]
        unique_adducts = tmp_table_1['pos_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['pos_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['neg_complexity'].sum() + tmp_table_2['pos_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['pos_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            neg_ion_list.append(tmp_table_1.loc[j, "neg_ion"])
            neg_adduct_list.append(tmp_table_1.loc[j, "neg_adduct"])
            neg_complexity_list.append(tmp_table_1.loc[j, "neg_complexity"])
            pos_ion_list.append(tmp_table_1.loc[j, "pos_ion"])
            pos_adduct_list.append(tmp_table_1.loc[j, "pos_adduct"])
            pos_complexity_list.append(tmp_table_1.loc[j, "pos_complexity"])
    candidates_table = list(zip(neg_ion_list, neg_adduct_list, neg_complexity_list, pos_ion_list, pos_adduct_list, pos_complexity_list))
    candidates_table = pd.DataFrame(candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity'])

    # Make a neutral table:
    neutral_idx = 0
    neutral_table = list()
    while len(candidates_table) > 0 :
        idx = candidates_table.index[0]
        neg_pool = list()
        pos_pool = list()
        idx_pool = list()
        neg_pool.append(candidates_table.loc[idx, "neg_ion"])
        pos_pool.append(candidates_table.loc[idx, "pos_ion"])
        new_neg_pool = neg_pool
        new_pos_pool = pos_pool
        idx_pool.append(idx)
        len_0 = 0
        while len_0 != len(idx_pool):
            len_0 = len(idx_pool)
            for i in new_neg_pool:
                idx_pool += candidates_table.index[candidates_table['neg_ion'] == i].tolist()
            for i in new_pos_pool:
                idx_pool += candidates_table.index[candidates_table['pos_ion'] == i].tolist()
            idx_pool = list(set(idx_pool))
        tmp_table_1 = candidates_table.loc[idx_pool]
        
        for i in tmp_table_1['neg_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['neg_ion'] == i][0]
            mgf_idx = int(merged_node_table.loc[i, "mgf_index"])
            ion_mz = merged_node_table.loc[i, "mz"]
            ion_rt = merged_node_table.loc[i, "rt"]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "neg_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_neg.index[adduct_table_base_neg["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_neg.loc[adduct_idx, "Adduct_code"]
            mol_mass = Neutral_mass_calculator(ion_mz,
                                adduct_table_base_neg.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_neg.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_neg.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "NEG")
            ion_rule_points = Species_rule(mgf_idx, mgf_file_neg)
            neutral_table.append((neutral_idx, mol_mass, ion_mz, ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "NEG"))
            
        for i in tmp_table_1['pos_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['pos_ion'] == i][0]
            mgf_idx = int(merged_node_table.loc[i, "mgf_index"])
            ion_mz = merged_node_table.loc[i, "mz"]
            ion_rt = merged_node_table.loc[i, "rt"]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "pos_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_pos.index[adduct_table_base_pos["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_pos.loc[adduct_idx, "Adduct_code"]
            mol_mass = Neutral_mass_calculator(ion_mz,
                                adduct_table_base_pos.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_pos.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_pos.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "POS")
            ion_rule_points = Species_rule(mgf_idx, mgf_file_pos)
            neutral_table.append((neutral_idx, mol_mass, ion_mz,  ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "POS"))            
        
        candidates_table.drop(idx_pool, inplace = True)
        neutral_idx += 1
    
    neutral_table = pd.DataFrame(neutral_table, columns = ['neutral_idx', 'neutral_mass', 'ion_mz', 'ion_rt', 'ion_idx', 'adduct', 'TIC', 'rule_points', 'samples', 'ion_mode'])         

    # Report the results: 
    print('Reporting results for NEG singletons paired to POS singletons...')
    for i in tqdm(neutral_table['neutral_idx'].unique()):
        tmp_table_1 = neutral_table[neutral_table['neutral_idx'] == i].copy()
        mol_mass = tmp_table_1['neutral_mass'].mean()
        mol_rt = tmp_table_1['ion_rt'].mean()
        mol_tic = tmp_table_1['TIC'].sum()
        mol_cluster = merged_node_table['cluster_id'].max() + 1
        adduct_count = len(tmp_table_1)
        tmp_table_1['mz_gap'] = abs(mol_mass - tmp_table_1['ion_mz'])
        tmp_table_1['rt_gap'] = abs(mol_rt - tmp_table_1['ion_rt'])
        mol_samples = '|'.join(tmp_table_1['samples'])
        mol_samples = list(set(mol_samples.split('|')))
        mol_samples.sort()
        mol_samples = '|'.join(mol_samples)
        new_idx = merged_node_table.index.max() + 1
        merged_node_table.loc[new_idx] = [mol_mass, mol_rt, mol_tic, 0, None, 0,
                             "neutral", adduct_count, None, "MIX", None, mol_samples, mol_cluster]
        for j in tmp_table_1.index:
            ion_idx = tmp_table_1.loc[j, 'ion_idx']
            adduct = tmp_table_1.loc[j, 'adduct']
            rule_points = tmp_table_1.loc[j, 'rule_points']
            mz_gap = tmp_table_1.loc[j, 'mz_gap']
            rt_gap = tmp_table_1.loc[j, 'rt_gap']
            tmp_mode = tmp_table_1.loc[j, "ion_mode"]
            merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
            merged_node_table.loc[ion_idx, 'status'] = "adduct"
            merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
            merged_node_table.loc[ion_idx, 'cluster_id'] = mol_cluster
            new_edge = merged_edge_table.index.max() + 1
            merged_edge_table.loc[new_edge] = [new_idx, ion_idx, 0 ,0 ,0 , rt_gap,
                                 mz_gap, "add_edge", None, adduct, adduct, tmp_mode]
            del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
            del_edge = del_edge[del_edge['node_2'] == ion_idx].index
            if len(del_edge) > 0 :
                merged_edge_table.drop(del_edge[0], inplace = True)
            


    # Connect singletons to singletons (only precursors and nodes from non molecular clusters)
    remains_table = cluster_table.index[~cluster_table["molecular_cluster"]].tolist()
    remains_table = [merged_node_table.index[merged_node_table['cluster_id'] == i].tolist() for i in remains_table]
    remains_table = list(flatten(remains_table)) # Added precursors and fragments from non-molecular clusters
    remains_table += merged_node_table.index[merged_node_table['cluster_id'] == -1].tolist() # Added singleton nodes
    precursor_ions = list(set(merged_node_table.index) - set(remains_table))
    precursor_ions = merged_node_table.loc[precursor_ions]
    precursor_ions = precursor_ions.index[precursor_ions['status'] == "precursor"].tolist()
    remains_table += precursor_ions # Added precursors from molecular clusters
    remains_table.sort()
    remains_table = merged_node_table.loc[remains_table]
    remains_table_pos = remains_table[remains_table['ion_mode'] == "POS"]
    remains_table_neg = remains_table[remains_table['ion_mode'] == "NEG"]

    print('Linking POS singletons to NEG singletons...')
    candidates_table = list()
    for i in tqdm(remains_table_pos.index):
        ion_1_mz = remains_table_pos.loc[i, "mz"]
        ion_1_rt = remains_table_pos.loc[i, "rt"]
        ion_1_samples = set(remains_table_pos.loc[i, "samples"].split('|'))
        neutral_table = Neutral_table(ion_1_mz, adduct_table_base_pos)
        hit_table_rt = remains_table_neg[remains_table_neg['rt'].between(ion_1_rt - rt_error, ion_1_rt + rt_error, inclusive = True)].copy()
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(ion_1_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        if len(hit_table_rt) == 0 : continue
        for j in neutral_table.index:
            ion_table = Ion_table(neutral_table.loc[j, "neutral_mass"], adduct_table_base_neg)
            for k in ion_table.index:
                hit_table = hit_table_rt[hit_table_rt['mz'].between(ion_table.loc[k, "ion_mz"] - mass_error, ion_table.loc[k, "ion_mz"] + mass_error, inclusive = True)]
                for l in hit_table.index:
                    candidates_table.append((i,
                                             neutral_table.loc[j, "Adduct"],
                                             neutral_table.loc[j, "Complexity"],
                                             l,
                                             ion_table.loc[k, "Adduct"],
                                             ion_table.loc[k, "Complexity"]))
    candidates_table = pd.DataFrame(candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity'])

    # Get the best hypotheses
    unique_pos = candidates_table['pos_ion'].unique().tolist()
    pos_ion_list = list()
    pos_adduct_list = list()
    pos_complexity_list = list()
    neg_ion_list = list()
    neg_adduct_list = list()
    neg_complexity_list = list()
    print('Selecting best POS annotations...')
    for i in tqdm(unique_pos) :
        tmp_table_1 = candidates_table[candidates_table['pos_ion'] == i]
        unique_adducts = tmp_table_1['pos_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['pos_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['pos_complexity'].sum() + tmp_table_2['neg_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['pos_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            pos_ion_list.append(tmp_table_1.loc[j, "pos_ion"])
            pos_adduct_list.append(tmp_table_1.loc[j, "pos_adduct"])
            pos_complexity_list.append(tmp_table_1.loc[j, "pos_complexity"])
            neg_ion_list.append(tmp_table_1.loc[j, "neg_ion"])
            neg_adduct_list.append(tmp_table_1.loc[j, "neg_adduct"])
            neg_complexity_list.append(tmp_table_1.loc[j, "neg_complexity"])
    candidates_table = list(zip(pos_ion_list, pos_adduct_list, pos_complexity_list, neg_ion_list, neg_adduct_list, neg_complexity_list))
    candidates_table = pd.DataFrame(candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity'])

    # Resolve NEG ions with multiple annotations
    unique_neg = candidates_table['neg_ion'].unique().tolist()
    pos_ion_list = list()
    pos_adduct_list = list()
    pos_complexity_list = list()
    neg_ion_list = list()
    neg_adduct_list = list()
    neg_complexity_list = list()
    print('Selecting best NEG annotations...')
    for i in tqdm(unique_neg):
        tmp_table_1 = candidates_table[candidates_table['neg_ion'] == i]
        unique_adducts = tmp_table_1['neg_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['neg_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['pos_complexity'].sum() + tmp_table_2['neg_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['neg_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            pos_ion_list.append(tmp_table_1.loc[j, "pos_ion"])
            pos_adduct_list.append(tmp_table_1.loc[j, "pos_adduct"])
            pos_complexity_list.append(tmp_table_1.loc[j, "pos_complexity"])
            neg_ion_list.append(tmp_table_1.loc[j, "neg_ion"])
            neg_adduct_list.append(tmp_table_1.loc[j, "neg_adduct"])
            neg_complexity_list.append(tmp_table_1.loc[j, "neg_complexity"])
    candidates_table = list(zip(pos_ion_list, pos_adduct_list, pos_complexity_list, neg_ion_list, neg_adduct_list, neg_complexity_list))
    candidates_table = pd.DataFrame(candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity'])

    # Make a neutral table:
    neutral_idx = 0
    neutral_table = list()
    while len(candidates_table) > 0 :
        idx = candidates_table.index[0]
        pos_pool = list()
        neg_pool = list()
        idx_pool = list()
        pos_pool.append(candidates_table.loc[idx, "pos_ion"])
        neg_pool.append(candidates_table.loc[idx, "neg_ion"])
        new_pos_pool = pos_pool
        new_neg_pool = neg_pool
        idx_pool.append(idx)
        len_0 = 0
        while len_0 != len(idx_pool):
            len_0 = len(idx_pool)
            for i in new_pos_pool:
                idx_pool += candidates_table.index[candidates_table['pos_ion'] == i].tolist()
            for i in new_neg_pool:
                idx_pool += candidates_table.index[candidates_table['neg_ion'] == i].tolist()
            idx_pool = list(set(idx_pool))
        tmp_table_1 = candidates_table.loc[idx_pool]
        
        for i in tmp_table_1['pos_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['pos_ion'] == i][0]
            mgf_idx = int(merged_node_table.loc[i, "mgf_index"])
            ion_mz = merged_node_table.loc[i, "mz"]
            ion_rt = merged_node_table.loc[i, "rt"]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "pos_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_pos.index[adduct_table_base_pos["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_pos.loc[adduct_idx, "Adduct_code"]
            mol_mass = Neutral_mass_calculator(ion_mz,
                                adduct_table_base_pos.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_pos.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_pos.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "POS")
            ion_rule_points = Species_rule(mgf_idx, mgf_file_pos)
            neutral_table.append((neutral_idx, mol_mass, ion_mz, ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "POS"))
            
        for i in tmp_table_1['neg_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['neg_ion'] == i][0]
            mgf_idx = int(merged_node_table.loc[i, "mgf_index"])
            ion_mz = merged_node_table.loc[i, "mz"]
            ion_rt = merged_node_table.loc[i, "rt"]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "neg_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_neg.index[adduct_table_base_neg["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_neg.loc[adduct_idx, "Adduct_code"]
            mol_mass = Neutral_mass_calculator(ion_mz,
                                adduct_table_base_neg.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_neg.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_neg.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "NEG")
            ion_rule_points = Species_rule(mgf_idx, mgf_file_neg)
            neutral_table.append((neutral_idx, mol_mass, ion_mz,  ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "NEG"))            
        
        candidates_table.drop(idx_pool, inplace = True)
        neutral_idx += 1
    
    neutral_table = pd.DataFrame(neutral_table, columns = ['neutral_idx', 'neutral_mass', 'ion_mz', 'ion_rt', 'ion_idx', 'adduct', 'TIC', 'rule_points', 'samples', 'ion_mode'])         

    # Report the results: 
    print('Reporting results for POS singletons paired to NEG singletons...')
    for i in tqdm(neutral_table['neutral_idx'].unique()):
        tmp_table_1 = neutral_table[neutral_table['neutral_idx'] == i].copy()
        mol_mass = tmp_table_1['neutral_mass'].mean()
        mol_rt = tmp_table_1['ion_rt'].mean()
        mol_tic = tmp_table_1['TIC'].sum()
        mol_cluster = merged_node_table['cluster_id'].max() + 1
        adduct_count = len(tmp_table_1)
        tmp_table_1['mz_gap'] = abs(mol_mass - tmp_table_1['ion_mz'])
        tmp_table_1['rt_gap'] = abs(mol_rt - tmp_table_1['ion_rt'])
        mol_samples = '|'.join(tmp_table_1['samples'])
        mol_samples = list(set(mol_samples.split('|')))
        mol_samples.sort()
        mol_samples = '|'.join(mol_samples)
        new_idx = merged_node_table.index.max() + 1
        merged_node_table.loc[new_idx] = [mol_mass, mol_rt, mol_tic, 0, None, 0,
                             "neutral", adduct_count, None, "MIX", None, mol_samples, mol_cluster]
        for j in tmp_table_1.index:
            ion_idx = tmp_table_1.loc[j, 'ion_idx']
            adduct = tmp_table_1.loc[j, 'adduct']
            rule_points = tmp_table_1.loc[j, 'rule_points']
            mz_gap = tmp_table_1.loc[j, 'mz_gap']
            rt_gap = tmp_table_1.loc[j, 'rt_gap']
            tmp_mode = tmp_table_1.loc[j, "ion_mode"]
            merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
            merged_node_table.loc[ion_idx, 'status'] = "adduct"
            merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
            merged_node_table.loc[ion_idx, 'cluster_id'] = mol_cluster
            new_edge = merged_edge_table.index.max() + 1
            merged_edge_table.loc[new_edge] = [new_idx, ion_idx, 0 ,0 ,0 , rt_gap,
                                 mz_gap, "add_edge", None, adduct, adduct, tmp_mode]
            del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
            del_edge = del_edge[del_edge['node_2'] == ion_idx].index
            if len(del_edge) > 0 :
                merged_edge_table.drop(del_edge[0], inplace = True)
    
    # Update status node nodes (pos neutrals, neg neutrals, pos adducts_ neg_adducts)
    merged_node_table.insert(merged_node_table.columns.get_loc('status') + 1, 
                             "status_universal", merged_node_table['status'].copy())
    merged_node_table['status'] = merged_node_table['ion_mode'].str.lower() + "_" + merged_node_table['status']
    
    # Update status for edges (pos add, neg add, pos frag, neg frag, pos single, neg single)
    merged_edge_table.insert(merged_edge_table.columns.get_loc('status') + 1, 
                             "status_universal", merged_edge_table['status'].copy())
    merged_edge_table['status'] = merged_edge_table['ion_mode'].str.lower() + "_" + merged_edge_table['status']
    
    # Round values:
    merged_node_table['mz'] = merged_node_table['mz'].round(4)
    merged_node_table['rt'] = merged_node_table['rt'].round(3)
    merged_edge_table['mz_gap'] = merged_edge_table['mz_gap'].round(4)
    merged_node_table.drop("samples", axis = 1, inplace = True)

    # Adding MZmine csv columns to the node table (sample intensities)   
    neutrals_idx = list(merged_node_table.index[merged_node_table['status'] == "neg_neutral"])
    neutrals_idx += list(merged_node_table.index[merged_node_table['status'] == "pos_neutral"])
    neutrals_idx += list(merged_node_table.index[merged_node_table['status'] == "mix_neutral"])
    pos_ions_idx = list(set(merged_node_table.index[merged_node_table['ion_mode'] == "POS"]) - set(neutrals_idx))
    neg_ions_idx = list(set(merged_node_table.index[merged_node_table['ion_mode'] == "NEG"]) - set(neutrals_idx))
    samples = list(neg_mzmine_csv.columns) + list(pos_mzmine_csv.columns)
    samples = list(set(samples))
    samples.sort()
    
    for sample in samples:
        merged_node_table[sample] = [0.0]*len(merged_node_table)
        
    print("Adding POS sample intensities to the merged node table...")
    for i in tqdm(pos_ions_idx):
        ion_id = merged_node_table.loc[i, "feature_id"]
        for sample in samples:
            merged_node_table.loc[i, sample] = pos_mzmine_csv.loc[ion_id, sample]
    
    print("Adding NEG sample intensities to the merged node table...")
    for i in tqdm(neg_ions_idx):
        ion_id = merged_node_table.loc[i, "feature_id"]
        for sample in samples:
            merged_node_table.loc[i, sample] = neg_mzmine_csv.loc[ion_id, sample]
    
    print("Adding neutral sample intensities to the merged node table...")
    for i in tqdm(neutrals_idx):
        ion_ids = list(merged_edge_table['node_2'][merged_edge_table['node_1'] == i])
        #ion_ids = merged_node_table.loc[ion_ids, "feature_id"]
        for sample in samples:
            merged_node_table.loc[i, sample] = merged_node_table.loc[ion_ids, sample].sum()
    
    
    #export the data
    if not os.path.isdir(out_full) :
        os.mkdir(out_full)
    if not os.path.isdir(out_samples) :
        os.mkdir(out_samples)
    
    merged_edge_table.to_csv(out_full + 'MIX_edges.csv', index_label = "Index")
    merged_node_table.to_csv(out_full + 'MIX_nodes.csv', index_label = "Index")
    
    if params['mm_export_samples'] : 
        Samplewise_export(neg_csv_file = in_path_full_neg + neg_csv,
                          pos_csv_file = in_path_full_pos + pos_csv,
                          out_path = out_samples,
                          merged_edge_table = merged_edge_table,
                          merged_node_table = merged_node_table)
    return



