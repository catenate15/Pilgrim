# level: HF
# Atomic number and non-scaled cartesian coordinates [bohr]
start_cc
   001   +0.00000000E+00  +0.00000000E+00  +0.00000000E+00
   006   +0.00000000E+00  +0.00000000E+00  +2.05169455E+00
   006   +2.74780919E+00  +0.00000000E+00  +3.04806777E+00
   001   -9.98174126E-01  +1.66790251E+00  +2.70788116E+00
   001   -1.01349036E+00  -1.66251490E+00  +2.69979691E+00
   008   +2.91878594E+00  +2.18641313E-01  +5.75050274E+00
   001   +3.73265129E+00  -1.69149385E+00  +2.36859028E+00
   001   +3.74759903E+00  +1.64217389E+00  +2.28627948E+00
   001   +2.12648802E+00  -1.50262707E+00  +6.43862761E+00
   001   +1.32624192E+00  -3.33489985E+00  +7.00707802E+00
end_cc

# Charge, multiplicity, energy [hartree],
# point group and rotational symmetry number
start_basic
   charge        0
   multiplicity  2
   energy       -152.58566390    # Total energy in hartree
   pointgroup    C1              # Point group
   rotsigma      1               # Rotational sigma
end_basic

# Non-scaled cartesian gradient [hartree/bohr]
start_grad
   +6.00000000E-08  +8.00000000E-08  +5.20000000E-07
   +2.00000000E-08  +4.31000000E-06  +1.15000000E-06
   -1.17000000E-06  -1.15000000E-06  +1.35000000E-06
   +2.19000000E-06  -3.55000000E-06  -1.40000000E-06
   -5.50000000E-07  -9.10000000E-07  +4.20000000E-07
   -8.80000000E-07  -1.72000000E-06  -1.67000000E-06
   -6.30000000E-07  +9.00000000E-07  +5.20000000E-07
   +2.20000000E-07  +4.50000000E-07  -9.00000000E-08
   +3.70000000E-07  +6.90000000E-07  -5.30000000E-07
   +3.80000000E-07  +9.10000000E-07  -2.90000000E-07
end_grad

# Low triangular part of force constant (hessian) matrix [hartree/bohr^2]
# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...
start_hess
   +7.09791500E-02  +3.56310000E-04  +6.88437500E-02  -1.98988000E-03  -2.10350000E-04
   +4.81364680E-01  -7.39443700E-02  +3.75200000E-05  +5.22773000E-03  +7.29586950E-01
   +2.80860000E-04  -7.96766700E-02  +4.78630000E-04  +2.16480000E-04  +8.61802310E-01
   +1.94111000E-03  +1.58080000E-04  -4.64109760E-01  -5.28963300E-02  +3.08744000E-03
   +8.45000650E-01  +4.41658000E-03  -4.93800000E-05  +3.81093000E-03  -2.91361330E-01
   +3.01070000E-04  -5.27096400E-02  +7.25107760E-01  +1.70400000E-05  -7.30700000E-04
   -2.50390000E-04  +2.28399000E-03  -1.27107130E-01  +3.14426000E-03  +1.98562000E-02
   +8.77295250E-01  -4.33516800E-02  +3.28000000E-04  -8.83490000E-03  -4.85675100E-02
   -7.74060000E-04  -1.25835450E-01  -4.77503100E-02  -1.58359000E-03  +7.61006340E-01
   +1.32323000E-03  -2.89991000E-03  -3.28857000E-03  -1.64768100E-01  +1.52073430E-01
   +6.22737700E-02  -1.81115000E-02  +3.60982300E-02  +1.60659000E-02  +1.67772920E-01
   -3.80559000E-03  +5.31192000E-03  +6.30616000E-03  +1.49612390E-01  -3.37528150E-01
   -1.01375910E-01  -4.08321000E-03  +5.98583000E-03  +2.80045000E-03  -1.61086960E-01
   +3.44222690E-01  +2.24606300E-02  -3.92225000E-02  -5.73111000E-03  +6.17173300E-02
   -1.00319610E-01  -1.19604050E-01  -9.43196000E-03  +1.65326800E-02  +6.71360000E-03
   -6.35793400E-02  +1.08507920E-01  +1.11815100E-01  +1.37184000E-03  +2.82479000E-03
   -3.24539000E-03  -1.67218100E-01  -1.52735730E-01  +6.21325300E-02  -1.79580700E-02
   -3.61487200E-02  +1.61050100E-02  +1.40519800E-02  +2.11653100E-02  -9.88202000E-03
   +1.70543940E-01  +3.61010000E-03  +4.93936000E-03  -6.30015000E-03  -1.50261430E-01
   -3.34038210E-01  +9.93584400E-02  +4.39164000E-03  +5.78705000E-03  -2.21106000E-03
   -2.16074100E-02  -1.96162500E-02  +1.53459600E-02  +1.62090540E-01  +3.41004960E-01
   +2.28106300E-02  +3.94822300E-02  -5.75445000E-03  +6.08447700E-02  +9.76486300E-02
   -1.18686230E-01  -8.94572000E-03  -1.61356400E-02  +6.35845000E-03  -9.83251000E-03
   -1.50142000E-02  +6.26226000E-03  -6.27730500E-02  -1.06942180E-01  +1.11749620E-01
   -8.43295000E-03  -4.29450000E-04  -1.76005000E-03  +3.87688000E-03  -3.83358000E-03
   -4.62017400E-02  -7.78143800E-02  -2.15412800E-02  +2.15918000E-03  +3.85843000E-03
   -4.35090000E-04  +1.72305000E-03  +1.93932000E-03  +1.95690000E-04  +7.54050000E-04
   +1.14244650E-01  -1.45031000E-03  +5.73130000E-04  +1.39410000E-04  -3.96790000E-04
   +5.22839000E-03  -4.65198000E-03  -1.79263600E-02  -1.17084140E-01  +1.93899400E-02
   -9.09940000E-04  -1.87070000E-04  -4.20900000E-05  +3.27980000E-04  -2.16160000E-04
   -2.78350000E-04  +1.08001270E-01  +2.95791370E-01  -3.79588000E-03  -8.25000000E-05
   +8.92380000E-04  -1.03144100E-02  -1.24470000E-04  -3.33282000E-02  -3.46255100E-02
   -6.74199200E-02  -3.97241900E-01  +1.80726000E-03  -2.63730000E-04  +1.36340000E-03
   +2.26000000E-03  +6.85000000E-05  +1.40234000E-03  +7.59039000E-03  -1.28154300E-02
   +5.49333250E-01  +1.86446000E-03  -4.91120000E-04  +6.08430000E-04  -1.82715200E-02
   +3.70014500E-02  +1.22321100E-02  -1.60027060E-01  +1.38935020E-01  +5.36587700E-02
   -6.54596000E-03  -1.71433000E-03  -3.70359000E-03  +2.92612000E-03  -7.40100000E-05
   +7.96280000E-04  +6.68028000E-03  -4.46750000E-03  +2.37710000E-04  +1.59831690E-01
   +1.93900000E-05  +1.40330000E-04  +2.58020000E-04  -5.69774000E-03  +5.33802000E-03
   +9.52450000E-04  +1.37575830E-01  -3.26407410E-01  -8.79940400E-02  -1.46598000E-03
   +8.75550000E-04  -6.84070000E-04  +4.97000000E-05  +9.01420000E-04  +1.51220000E-04
   -1.12071000E-03  +3.22843000E-03  -3.10191000E-03  -1.49130160E-01  +3.32516510E-01
   +1.19865000E-03  +5.30400000E-05  +9.16990000E-04  -8.18337000E-03  +1.40372600E-02
   +8.07991000E-03  +5.79176100E-02  -8.97370100E-02  -1.05577820E-01  -3.78089000E-03
   -1.02148000E-03  -8.39670000E-04  +2.06940000E-04  -3.11540000E-04  -3.33600000E-05
   +1.61446800E-02  -4.24512000E-02  -3.24018200E-02  -5.53438400E-02  +1.02823500E-01
   +1.20198270E-01  +2.18472000E-03  +4.51380000E-04  +6.17120000E-04  -1.84716800E-02
   -3.52465400E-02  +1.41289600E-02  -1.61863000E-01  -1.41326400E-01  +6.57255800E-02
   +2.70207000E-03  +1.34940000E-04  +7.75550000E-04  -6.33195000E-03  +1.55444000E-03
   -3.39996000E-03  +3.82100000E-03  +2.69502000E-03  +1.02813000E-03  +1.33420700E-02
   +2.04535700E-02  -9.25894000E-03  +1.63226670E-01  +2.98600000E-05  +2.13490000E-04
   -1.78410000E-04  +5.03737000E-03  +5.47206000E-03  -1.54489000E-03  -1.38816190E-01
   -3.18982230E-01  +1.02398970E-01  -5.63000000E-05  +8.64660000E-04  -1.70270000E-04
   +1.41118000E-03  +9.07880000E-04  +6.11350000E-04  +2.89223000E-03  +7.02675000E-03
   +2.70837000E-03  -2.11460400E-02  -1.97853600E-02  +1.45226000E-02  +1.50158620E-01
   +3.21807120E-01  +7.87620000E-04  -2.46460000E-04  +9.30730000E-04  -8.94366000E-03
   -1.40022800E-02  +7.26859000E-03  +6.50697900E-02  +1.02128250E-01  -1.19103170E-01
   +3.49500000E-04  +2.55000000E-04  +4.56200000E-05  -3.84175000E-03  +9.45510000E-04
   -7.83630000E-04  +2.06528700E-02  +3.85298300E-02  -1.97628000E-02  -7.87969000E-03
   -1.26700300E-02  +7.70541000E-03  -6.72915600E-02  -1.15686740E-01  +1.33510360E-01
   +5.63850000E-04  +2.31490000E-04  +6.39600000E-05  +2.74250000E-04  +1.91488000E-03
   -4.58600000E-05  -2.66990000E-03  +1.31150000E-03  -1.77909000E-02  -3.21320000E-04
   +2.64730000E-04  -1.06320000E-04  +8.45110000E-04  +1.30900000E-05  +6.90000000E-06
   -2.84065600E-02  -2.98885700E-02  +1.30044300E-02  +3.92900000E-04  -1.77740000E-04
   +1.96038000E-03  +1.78309000E-03  +6.86960000E-04  +1.83693000E-03  +1.04434700E-02
   +1.71259000E-03  +6.38540000E-04  -2.75190000E-04  -1.71935000E-03  +9.17720000E-04
   +1.78024000E-03  -3.63998000E-03  -1.29200000E-04  -4.02060300E-02  -2.14710000E-04
   +1.85370000E-04  +6.25200000E-05  +1.54523000E-03  +5.51790000E-04  +4.34580000E-04
   -2.51383800E-02  -7.18533400E-02  +2.87791800E-02  +1.53593000E-03  +3.79570000E-03
   +4.21314000E-03  +1.98440000E-03  +3.14225000E-03  +2.48755000E-03  -3.82390400E-02
   -6.99638800E-02  -3.85210000E-04  -4.47930000E-04  +3.54910000E-04  +1.52183000E-03
   -2.75310000E-04  +1.07318000E-03  +2.68091000E-02  +5.28558800E-02  -1.44881300E-02
   +2.84200000E-05  -2.27470000E-04  -3.50700000E-05  -1.13848000E-03  -2.95910000E-04
   -2.29130000E-04  -1.90971200E-02  -3.76777200E-02  -5.80739000E-02  -6.61730000E-04
   -1.81180000E-04  +1.15549000E-03  -2.51232000E-03  -2.74901000E-03  -1.05325100E-02
   +2.21543000E-02  +5.04226900E-02  +7.14968800E-02  -3.26520000E-04  -3.16400000E-05
   -4.42700000E-05  +2.97030000E-04  +2.76900000E-05  -8.54910000E-04  +2.80910000E-04
   +5.14410000E-04  +3.74596000E-03  +3.82600000E-05  -5.21900000E-05  +2.66700000E-05
   -1.70190000E-04  +8.73400000E-05  -2.61390000E-04  -1.97666700E-02  -5.59848000E-02
   +2.28078700E-02  -1.92980000E-04  -5.06160000E-04  -8.61220000E-04  -3.93000000E-04
   -1.97700000E-04  -7.40060000E-04  +1.70951200E-02  +6.21733000E-02  -2.67188000E-02
   +3.13804000E-03  -7.70250000E-04  -2.53150000E-04  +3.22800000E-05  +8.87550000E-04
   -4.08330000E-04  -9.08130000E-04  +2.39038000E-03  +1.37267000E-03  +7.85141000E-03
   +6.95500000E-05  -1.14540000E-04  -1.05300000E-05  -5.30290000E-04  -2.21840000E-04
   +4.23600000E-05  -5.85907100E-02  -1.22507370E-01  +5.22519100E-02  -4.49240000E-04
   -6.03190000E-04  -2.12832000E-03  -8.59440000E-04  -6.66620000E-04  -1.74063000E-03
   +6.38826900E-02  +1.32715060E-01  -6.14240500E-02  -6.03025000E-03  -9.31270000E-03
   +3.24010000E-04  +1.88400000E-04  -2.94600000E-05  -4.06390000E-04  +2.43780000E-04
   +1.41360000E-04  -1.44280000E-04  +4.65470000E-04  -2.99703000E-03  -4.35500000E-05
   +3.32500000E-05  +9.91000000E-06  +1.76200000E-04  +3.42430000E-04  -2.85870000E-04
   +1.80346900E-02  +3.98575800E-02  -1.21827600E-02  +5.55500000E-05  +4.46030000E-04
   +7.96590000E-04  +1.87430000E-04  +8.80300000E-05  +7.21400000E-04  -2.10838300E-02
   -4.76986800E-02  +9.27826000E-03  +2.90016000E-03  +6.03370000E-03  +4.54759000E-03
end_hess
