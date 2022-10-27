ffligandsString = ['gaff', 'gaff2', 'cgenff']

standardChargeMethods = ['am1bcc', 'am1bcc-pol', 'mmff94', 'eem', 'qeq', 'qtpie']
standardChargeMethodsHelp = ['am1bcc: AM1 with bond charge correction (antechamber)',
                     'am1bcc-pol: STaGE\'s own more polarized bond charge correction (antechamber)',
                     'mmff94: MMFF94 (Open Babel)',
                     'eem: electronegativity equalization method (Open Babel)',
                     'qeq: Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991) (Open Babel)',
                     'qtpie: Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007) (Open Babel)']
extraChargeMethods = ['gaussian/hf']
extraChargeMethodsHelp = ['gaussian/hf: Hatree-Fock/6-31G(d) basis set followed by RESP (Gaussian)']
                          #'gamess/hf: Hatree-Fock/6-31G(d) basis set followed by RESP (GAMESS)',
chargeMethods = standardChargeMethods + extraChargeMethods
chargeMethodsHelp = standardChargeMethodsHelp + extraChargeMethodsHelp

waterModels = ['opc', 'spce', 'tip4pew', 'spc', 'tip3p']

boxTypes = ['triclinic', 'cubic', 'dodecahedron', 'octahedron']
