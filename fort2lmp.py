#  =========================================================
#  Convert Fortran generated topology to LAMMPS input format
#                      September, 2020
#  =========================================================

# input params.
aa = 33.31
bb = 33.31
cc = bb
sigma = 3.4

# read F90 file
o = open('topF90.Ar', 'r')
crd = []

for line in o:
    atm = []
    b = line.split(' ')
    for elt in b:
        if (elt != ''):
            atm.append(elt)
    crd.append(atm)

# create LAMMPS .data format
o = open('top.Ar', 'w')
nUnit = len(crd)
o.write('Argon random initial configuration'+'\n')
o.write('\n')
o.write(str(nUnit) + ' ' + 'atoms' + '\n')
o.write('0 bonds' + '\n')
o.write('0 angles' + '\n')
o.write('0 dihedrals' + '\n')
o.write('0 impropers' + '\n')
o.write('\n')
o.write('1 atom types' + '\n')
o.write('\n')
o.write('0.0 ' + str(aa/sigma) + ' xlo xhi'+'\n')
o.write('0.0 ' + str(bb/sigma) + ' ylo yhi'+'\n')
o.write('0.0 ' + str(cc/sigma) + ' zlo zhi'+'\n')
o.write('\n')
o.write('Masses' + '\n')
o.write('\n')
o.write('  1 ' + str(1.0) + '\n')
o.write('\n')
o.write('Atoms' + '\n')
o.write('\n')

for at in crd:
    x = float(at[1])/sigma
    y = float(at[2])/sigma
    z = float(at[3])/sigma
    o.write('  ' + str(at[0]) + ' 1  ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n')

o.close()
