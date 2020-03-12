# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:54:07 2020

@author: gawe
"""

# refractive focusing elements - Paul F. Goldsmith - Quasioptical Systems pg 81

def _mats(material, metal=False):
    # ============== Vacuum filled guide ============== #
    if material.upper().find("VACUUM")>-1 or material.upper().find("NONE")>-1 or material is None:
        epsr = 1.0
        loss_tangent = 0.0

    # ============== Air filled guide ============== #
    elif material.upper().find("AIR")>-1:
        epsr = 1.0006  # relative permeability of material in the guide, air
        loss_tangent = 0.0  # loss tangent of material in guide

    # ============== Macor ceramic ============== #
    elif material.lower().find("macor")>-1:
#        nindex = 2.38
#        epsr = nindex**2.0
#        loss_tangent = 275e-4

#        fref = lower
        epsr = 5.67  # relative permeability of material in the guide
        loss_tangent = 0.0071  # loss tangent of material in guide

    # ============== Magnesium Oxide ============== #
    elif material.lower().find("magnesium")>-1:
#        fref = 92.8e9
        nindex = 3.132
        epsr = nindex**2.0
        loss_tangent = 0.46e-4

    # ============== Mica filled guide ============== #
    elif material.lower().find("mica")>-1 or material.lower().find("glimmer")>-1:
#        fref = 120-1000e9
        nindex = 2.54
        epsr = nindex**2.0
        loss_tangent = 13e-4
#        epsr = 5.7  # relative permeability of material in the guide
#        loss_tangent = 0.000  # loss tangent of material in guide

    # ============== Mylar ============== #
    elif material.lower().find("mylar")>-1:
#        fref = 120-1000e9
        nindex = 1.73
        epsr = nindex**2.0
        loss_tangent = 360e-4

#        fref = 140e9
#        nindex = 1.83
#        epsr = nindex**2.0
#        loss_tangent = 100 + 100/-50e-4

    # ============== Paraffin ============== #
    elif material.lower().find("paraffin")>-1:
#        fref = 289e9
        nindex = 1.50
        epsr = nindex**2.0
        loss_tangent = 34e-4

    # ============== Nylon ============== #
    elif material.lower().find("nylon")>-1:
#        fref = 130-180e9
        nindex = 1.7267
        epsr = nindex**2.0
        loss_tangent = 96e-4

    # ============== Plexiglas ============== #
    elif material.lower().find("plex")>-1:
#        fref = 150-6000e9
        nindex = 1.6067
        epsr = nindex**2.0
        loss_tangent = 87e-4

    # ============== Polyethylene ============== #
    elif material.lower().find("polyeth")>-1 or material.lower().find("pe")>-1:
#        fref = 90-270e9
        nindex = 1.519
        epsr = nindex**2.0
        loss_tangent = 3.6e-4

    # ============== Polypropylene ============== #
    elif material.lower().find("polypro")>-1 or material.lower().find("pp")>-1:
#        fref = 90-270e9
        nindex = 1.502
        epsr = nindex**2.0
        loss_tangent = 5.6e-4

    # ============== Polystyrene filled guide ============== #
    elif material.lower().find("polystyrene")>-1 or material.lower().find('ps')>-1:
#        fref = 140e9
        nindex = 1.59
        epsr = nindex**2.0
        loss_tangent = 20e-4

#        epsr = 2.4  # relative permeability of material in the guide, air
#        loss_tangent = 0.0  # loss tangent of material in guide

    # ============== Pyrex ============== #
    elif material.lower().find("pyrex")>-1:
#        fref = 250-400e9
        nindex = 2.11
        epsr = nindex**2.0
        loss_tangent = 28e-4

    # ============== Fused Quartz filled guide ============== #
    elif material.lower().find("quartz")>-1:
#        fref = 60e9
        nindex = 2.106
        epsr = nindex**2.0
        loss_tangent = 0.6e-4

#        epsr = 3.78  # relative permeability of material in the guide
#        loss_tangent = 0.000  # loss tangent of material in guide

    # ============== Rexolite ============== #
    elif material.lower().find("rexolite")>-1:
#        fref = 140e9
        nindex = 1.57
        epsr = nindex**2.0
        loss_tangent = 20e-4

    # ============== Sapphire filled guide ============== #
    elif material.lower().find("saphire")>-1:
#        fref = 180e9
        nindex = 3.094
        epsr = nindex**2.0
        loss_tangent = 5.8e-4

#        epsr = 10.0  # relative permeability of material in the guide
#        loss_tangent = 0.000  # loss tangent of material in guide

    # ============== Silicon filled guide ============== #
    elif material.lower().find("silicon")>-1 or material.lower().find("si")>-1:
#        fref = 90-450e9
        nindex = 3.417
        epsr = nindex**2.0
        loss_tangent = 6e-4

    # ============== Spectralon ============== #
    elif material.lower().find("spectralon")>-1:
#        fref = 291e9
        nindex = 1.31
        epsr = nindex**2.0
        loss_tangent = 213e-4

    # ============== Spinel ============== #
    elif material.lower().find("spinel")>-1:
#        fref = 90-350e9
        nindex = 2.8942
        epsr = nindex**2.0
        loss_tangent = 5e-4

    # ============== Styrofoam ============== #
    elif material.lower().find("styro")>-1:
#        fref = 245e9
        nindex = 1.017
        epsr = nindex**2.0
        loss_tangent = 0.53e-4

    # ============== Teflon (PTFE) filled guide ============== #
    elif material.lower().find("ptfe")>-1 or material.lower().find("teflon")>-1:
#        fref = 90-270e9
        nindex = 1.43855
        epsr = nindex**2.0
        loss_tangent = 5.3e-4

#        epsr = 2.1  # relative permeability of material in the guide
#        loss_tangent = 0.001  # loss tangent of material in guide

    # ============== Titanium Dioxide ============== #
    elif (material.lower().find("titan")>-1 and material.lower().find("diox")>-1) or material.lower().find("tidiox")>-1:
#        fref = 10.125e9
        nindex = 9.54
        epsr = nindex**2.0
        loss_tangent = 5e-4

    # ============== Polymethylpentene ============== #
    elif material.lower().find("polymethyl")>-1 or material.lower().find("tpx")>-1:  # is this actually for 'polyimide'?
#        fref = 70-270e9
        nindex = 1.458
        epsr = nindex**2.0
        loss_tangent = 5e-4

#        epsr = 4.3  # relative permeability of material in the guide
#        loss_tangent = 0.004  # loss tangent of material in guide

    # ============== Zinc Sulfide ============== #
    elif (material.lower().find("zinc")>-1 and material.lower().find("sulf")>-1) or material.lower().find("zns")>-1:
#        fref = 36e9
        nindex = 2.93
        epsr = nindex**2.0
        loss_tangent = 13e-4

    # ============== Zinc Selenide ============== #
    elif (material.lower().find("zinc")>-1 and material.lower().find("sele")>-1) or material.lower().find("znse")>-1:
#        fref = 891e9
        nindex = 3.125
        epsr = nindex**2.0
        loss_tangent = 33.1e-4

    # ============== Polyamide filled guide ============== #
    elif material.lower().find("polyamide")>-1 or material.lower().find("kapton")>-1:  # is this actually for 'polyimide'?
        epsr = 4.3  # relative permeability of material in the guide
        loss_tangent = 0.004  # loss tangent of material in guide

    # ============== Alumina ceramic ============== #
    elif material.lower().find("alumina")>-1:
        epsr = 8.66  # relative permeability of material in the guide
        loss_tangent = 0.0018  # loss tangent of material in guide

    # endif

    if metal:
        if material.lower().find("copper")>-1 or material.lower().find('cu')>-1:
            mur = 0.999991 # Copper, relative permeability of cavity walls
            rho = 1.724e-8 # Copper, ohm-meter, resistivity of walls of cavity
        elif material.lower().find("aluminum")>-1 or material.lower().find('al')>-1:
            mur = 1.00002 # Aluminum, relative permeability of cavity walls
            rho = 2.65e-8 # Aluminum, ohm-meter, resistivity of walls of cavity
        elif material.lower().find("brass")>-1 or material.lower().find('bra')>-1:
            mur = 1.05 # Brass, relative permeability of cavity walls
            rho = 6.39e-8 # Brass, ohm-meter, resistivity of walls of cavity
        # end if
        return mur, rho
    # end if
    return epsr, loss_tangent
# end def materials
